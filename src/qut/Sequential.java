package qut;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class Sequential {
    private static ThreadLocal<Series> threadLocalPattern = ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static ConcurrentHashMap<String, Sigma70Consensus> consensus = new ConcurrentHashMap<>();
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

    private static int THREAD_COUNT = Runtime.getRuntime().availableProcessors();
    private static ExecutorService threadPool = Executors.newFixedThreadPool(THREAD_COUNT);

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
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location - 1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
        else {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i = 0; i < upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart - i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(threadLocalPattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir) {
        // Commented out code is used to allow the program to work with 4 or less threads
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

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException {

        // Start time for timer
        long startTime = System.currentTimeMillis();
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<String> fileNames = ListGenbankFiles(dir);

        List<Future<?>> futures = new ArrayList<>();
        for (String filename : fileNames) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                futures.add(threadPool.submit(() -> {
                    System.out.println(referenceGene.name);
                    for (Gene gene : record.genes)
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
                }));
            }
        }

        // Wait for all tasks to complete
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        // Shutdown the thread pool
        threadPool.shutdown();

        // End time for timer of code
        long endTime = System.currentTimeMillis();

        // Get time taken to complete program and display it
        long duration = endTime - startTime;
        System.out.println("Time taken: " + duration + " milliseconds");

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        System.out.println(THREAD_COUNT);
        // Run the program
        run("./referenceGenes.list", "./Ecoli");
    }
}
