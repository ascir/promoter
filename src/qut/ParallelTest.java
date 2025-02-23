package qut;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class ParallelTest {

    private Parallel parallel;

    @BeforeEach
    public void setUp() {
        parallel = new Parallel();
    }

    @Test
    public void testParallelStreams() throws IOException {
        // Declare default consensus results
        HashMap<String, String> sequentialResult = new HashMap<>();
        sequentialResult.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        sequentialResult.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        sequentialResult.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        sequentialResult.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        sequentialResult.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        sequentialResult.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        sequentialResult.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        sequentialResult.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        sequentialResult.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");

        // Run the original Sequential version
        parallel.run("../referenceGenes.list", "../Ecoli");

        // Run the Parallel version with the same input dataset
        parallel.runParallelStreams("../referenceGenes.list", "../Ecoli", 8);

        // Get the consensus results
        HashMap<String, Sigma70Consensus> parallelConsensus = parallel.getConsensus();

        // Compare the consensus results
        compareConsensusResults(sequentialResult, parallelConsensus);
    }

    @Test
    public void testExplicitThreads() throws IOException {
        // Declare default consensus results
        HashMap<String, String> sequentialResult = new HashMap<>();
        sequentialResult.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        sequentialResult.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        sequentialResult.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        sequentialResult.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        sequentialResult.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        sequentialResult.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        sequentialResult.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        sequentialResult.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        sequentialResult.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");

        // Run the Parallel version with explicit threads
        parallel.runExplicitThreads("../referenceGenes.list", "../Ecoli");

        // Get the consensus results
        HashMap<String, Sigma70Consensus> parallelConsensus = parallel.getConsensus();

        // Compare the consensus results
        compareConsensusResults(sequentialResult, parallelConsensus);
    }

    @Test
    public void testExecutorService() throws IOException, ExecutionException, InterruptedException {
        // Declare default consensus results
        HashMap<String, String> sequentialResult = new HashMap<>();
        sequentialResult.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        sequentialResult.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        sequentialResult.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        sequentialResult.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        sequentialResult.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        sequentialResult.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        sequentialResult.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        sequentialResult.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        sequentialResult.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");

        // Run the Parallel version with ExecutorService
        parallel.runExecutorService("../referenceGenes.list", "../Ecoli", 8);
    
        // Get the consensus results
        HashMap<String, Sigma70Consensus> parallelConsensus = parallel.getConsensus();

        // Compare the consensus results
        compareConsensusResults(sequentialResult, parallelConsensus);
    }

    // Compare the consensus results with the default results
    public void compareConsensusResults(HashMap<String, String> sequentialResult,
                                       HashMap<String, Sigma70Consensus> parallelConsensus) {
        // Compare each consensus result
        for (Map.Entry<String, String> entry : sequentialResult.entrySet()) {
            String geneName = entry.getKey();
            String defaultResult = entry.getValue();
            String parallelResult = parallelConsensus.get(geneName).toString();

            assertEquals(defaultResult, parallelResult, "Consensus result mismatch for " + geneName);
        }
    }
}
