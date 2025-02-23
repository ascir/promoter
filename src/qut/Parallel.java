package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class Parallel {
	private static HashMap<String, Sigma70Consensus> consensus = new HashMap<>();
//        private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
	private static final ThreadLocal<Series> sigma70_pattern = ThreadLocal
			.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
	private static final Matrix BLOSUM_62 = BLOSUM62.Load();
	private static byte[] complement = new byte['z'];

	static {
		complement['C'] = 'G';
		complement['c'] = 'g';
		complement['G'] = 'C';
		complement['g'] = 'c';
		complement['T'] = 'A';
		complement['t'] = 'a';
		complement['A'] = 'T';
		complement['a'] = 't';
	}

	public static HashMap<String, Sigma70Consensus> getConsensus() {
		return consensus;
	}

	private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
		List<Gene> referenceGenes = new ArrayList<Gene>();
		while (true) {
			String name = reader.readLine();
			if (name == null)
				break;
			String sequence = reader.readLine();
			referenceGenes.add(new Gene(name, 0, 0, sequence));
			consensus.put(name, new Sigma70Consensus());
		}
		consensus.put("all", new Sigma70Consensus());
		reader.close();
		return referenceGenes;
	}

	private static boolean Homologous(PeptideSequence A, PeptideSequence B) {
		return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f)
				.calculateScore() >= 60;
	}

	private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
		int upStreamDistance = 250;
		if (gene.location < upStreamDistance)
			upStreamDistance = gene.location - 1;

		if (gene.strand == 1)
			return new NucleotideSequence(
					Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
		else {
			byte[] result = new byte[upStreamDistance];
			int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
			for (int i = 0; i < upStreamDistance; i++)
				result[i] = complement[dna.bytes[reverseStart - i]];
			return new NucleotideSequence(result);
		}
	}

	private static Match PredictPromoter(NucleotideSequence upStreamRegion) {
		return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
	}

	private static void ProcessDir(List<String> list, File dir) {
		if (dir.exists())
			for (File file : dir.listFiles())
				if (file.isDirectory())
					ProcessDir(list, file);
				else
					list.add(file.getPath());
	}

	private static List<String> ListGenbankFiles(String dir) {
		List<String> list = new ArrayList<String>();
		ProcessDir(list, new File(dir));
		return list;
	}

	private static GenbankRecord Parse(String file) throws IOException {
		GenbankRecord record = new GenbankRecord();
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		record.Parse(reader);
		reader.close();
		return record;
	}

	public synchronized void addConsensus(String name, Match prediction) {
		consensus.get(name).addMatch(prediction);
		consensus.get("all").addMatch(prediction);
	}

	/*
	 * =============================================================================
	 * 						 PARALLEL WITH EXPLICIT THREADS
	 * =============================================================================
	 * 
	 */

    private static class FileProcessorTask implements Runnable {
        private String filename;
        private List<Gene> referenceGenes;

        public FileProcessorTask(String filename, List<Gene> referenceGenes) {
            this.filename = filename;
            this.referenceGenes = referenceGenes;
        }

        @Override
        public void run() {
            try {
                System.out.println(filename);
                GenbankRecord record = Parse(filename);
                for (Gene referenceGene : referenceGenes) {
                    System.out.println(referenceGene.name);
                    for (Gene gene : record.genes) {
                        if (Homologous(gene.sequence, referenceGene.sequence)) {
                            NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                            Match prediction = PredictPromoter(upStreamRegion);
                            if (prediction != null) {
                                synchronized (consensus) {
                                    consensus.get(referenceGene.name).addMatch(prediction);
                                    consensus.get("all").addMatch(prediction);
                                }
                            }
                        }
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public static void runExplicitThreads(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<String> files = ListGenbankFiles(dir);
        List<Thread> threads = new ArrayList<>();

        for (String filename : files) {
            Thread thread = new Thread(new FileProcessorTask(filename, referenceGenes));
            thread.start();
            threads.add(thread);
        }

        // Wait for all threads to complete
        for (Thread thread : threads) {
            try {
                thread.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // Output the consensus results
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }


	/*
	 * **===========================================================================
	 * 									 PARALLEL STREAMS
	 ** ============================================================================
	 */

	public void runParallelStreams(String referenceFile, String dir, int numOfThreads) throws IOException {
		// Initialize a list to hold FileProcessor tasks
		List<FileProcessor> fprocessor = new ArrayList<>();
		// Parse reference genes from the reference file
		List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
		// Iterate over each reference gene and GenBank file
		for (Gene referenceGene : referenceGenes) {
			for (String filename : ListGenbankFiles(dir)) {
				GenbankRecord record = Parse(filename);
				for (Gene gene : record.genes) {
					fprocessor.add(new FileProcessor(referenceGene, gene, record));
				}
			}
		}
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(numOfThreads));
		fprocessor.parallelStream()
				.filter(task -> Homologous(task.getGene().sequence, task.getReferenceGene().sequence)).forEach(task -> {
					NucleotideSequence upStreamRegion = GetUpstreamRegion(task.getRecord().nucleotides, task.getGene());
					Match prediction = PredictPromoter(upStreamRegion);
					if (prediction != null)
						addConsensus(task.getReferenceGene().name, prediction);
				});

		for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
			System.out.println(entry.getKey() + " " + entry.getValue());
	}

	/*
	 * =============================================================================
	 							EXECUTOR SERVICE CODE
	 * =============================================================================
	 */
	
	public void runExecutorService(String referenceFile, String dir, int threadNum)
			throws IOException, ExecutionException, InterruptedException {
		// ExecutorService executorService =
		// Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		ExecutorService executorService = Executors.newFixedThreadPool(threadNum);
		List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
		List<Future> futureTasks = new ArrayList<>();

		for (Gene referenceGene : referenceGenes) {
			for (String filename : ListGenbankFiles(dir)) {
				GenbankRecord record = Parse(filename);
				for (Gene gene : record.genes) {
					Future futureTask = executorService.submit(new RunnableTask(referenceGene, gene, record));
					futureTasks.add(futureTask);
				}
			}
		}


		executorService.shutdown();
		for (Future futureTask : futureTasks)
			futureTask.get();
		for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
			System.out.println(entry.getKey() + " " + entry.getValue());
	}

	public class RunnableTask implements Runnable {
		private final Gene referenceGene;
		private final Gene gene;
		private final GenbankRecord record;

		public RunnableTask(Gene referenceGene, Gene gene, GenbankRecord record) {
			this.referenceGene = referenceGene;
			this.gene = gene;
			this.record = record;
		}

		@Override
		public void run() {
			if (Homologous(gene.sequence, referenceGene.sequence)) {
				NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
				Match prediction = PredictPromoter(upStreamRegion);
				if (prediction != null) {
					addConsensus(referenceGene.name, prediction);
				}
			}
		}
	}

	/*
	 * =============================================================================
	  									SEQUENTIAL
	 * =============================================================================
	 */

	public static void run(String referenceFile, String dir) throws IOException {
		List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
		// Loop Transformation: for loop 1 <--> for loop 2
		for (Gene referenceGene : referenceGenes) {
			System.out.println(referenceGene.name);
			for (String filename : ListGenbankFiles(dir)) {
				System.out.println(filename);
				GenbankRecord record = null;
				try {
					record = Parse(filename);
				} catch (IOException e) {
					e.printStackTrace();
				}
				for (Gene gene : record.genes) {
					if (Homologous(gene.sequence, referenceGene.sequence)) {
						NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
						Match prediction = PredictPromoter(upStreamRegion);
						if (prediction != null) {
							consensus.get(referenceGene.name).addMatch(prediction);
							consensus.get("all").addMatch(prediction);
						}
					}
				}
			}
		}
	}

	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		long startTime = System.currentTimeMillis();

		// Uncomment 1 method you want to run
		// Explicit Threads
//        runExplicitThreads("../referenceGenes.list", "../Ecoli");
		// Parallel Stream
//		new Parallel().runParallelStreams("../referenceGenes.list", "../Ecoli", 8);
		// Executor Service
        new Parallel().runExecutorService("../referenceGenes.list", "../Ecoli", 8);
		// Sequential run
//        run("src/referenceGenes.list", "src/Ecoli");

		long execTime = System.currentTimeMillis() - startTime;
		System.out.println("Execution time: " + execTime / 1000 + " s");
	}
}