This aim of this application is to identify and analyse promoter regions in bacterial genes by comparing them with a set of reference genes. This includes steps such as parsing reference genes and their sequences, and check if it is homologous to any gene in the genbank files using the Smith-Waterman-Gotoh algorithm. If homologous, the promotor region in the upstream DNA is predicted using the Sigma70 consensus sequence and a consensus map is updated and displayed.

### Process

•	Parse reference genes and their sequences.
•	Parse GenBank files containing E. coli gene sequences.
•	For each gene in the GenBank files, check if it's homologous (shares common evolutionary origin) to any reference gene.
•	If homologous, predict the promoter region in the upstream DNA sequence of the gene.
•	Update a consensus map with the prediction.

### Describing the codebase

The Sequential.java class is the main entry point for the application, which hosts the main method. The application contains the following methods for different functionalities:
•	Main: calls a run method with a reference gene list and a directory containing E.coli GenBank files. 
•	run: This is the core method that processes reference genes and GenBank files. This method contains 4 for loops, of which 3 are nested. 
•	The first, outermost for loop iterates over the GenBank files, which are quite extensive (up to 20mb). For each record, it parses the file to extract genomic data into a GenbankRecord object.
•	The second loop iterates over the list of reference genes and checks for homology with the current GenBank record.
•	The third, innermost loop iterates over all the genes in the record and checks if the gene is homologous to the current reference gene using the Homologous method. If found to be homologous, it predicts the promoter region for that gene and updates the consensus map with the prediction.
•	The fourth independent for loop iterates over the consensus map and prints the name of each reference gene and its associated Sigma70 consensus prediction.
•	Parse: Parses a GenBank file and returns a GenbankRecord object.
•	ParseReferenceGenes: Reads a reference file and returns a list of reference genes. It also populates the consensus map with Sigma70Consensus objects for each gene.
•	Homologous: Determines if two peptide sequences are homologous using the Smith-Waterman-Gotoh alignment algorithm.
•	GetUpstreamRegion: Retrieves the upstream region of a gene from a given DNA sequence.
•	PredictPromoter: Predicts the promoter region in an upstream DNA sequence using the Sigma70 pattern.
•	ProcessDir & ListGenbankFiles: Utility methods to list all GenBank files in a given directory.
