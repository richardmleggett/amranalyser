/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package amranalyser;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * @author leggettr
 */
public class AMRAnalyser {
    private ArrayList<BlastSet> samples = new ArrayList<BlastSet>();
    private ArrayList<GeneSet> geneSets = new ArrayList<GeneSet>();
    private ArrayList<String> sampleIds = new ArrayList<String>();
    private CardGeneParser cgp = new CardGeneParser("/Users/leggettr/Desktop/BAMBI/all_genes_grouping_70pc.txt");

    public AMRAnalyser() {
    }

    private void processNewSample(String sampleId, BlastSet sample) {
        GeneSet gs = new GeneSet(sampleId);

        cgp.parseBlastAlignments(sample, gs, null);

        sampleIds.add(sampleId);
        geneSets.add(gs);
        samples.add(sample);
    }
    
    public void addSampleFromChunkSet(String sampleId, double e, int l, int i, String directory, String prefix, String midfix) {
        BlastSet sample = new BlastSet(e, l, i);
        sample.parseChunkSet(directory, prefix, midfix);
        processNewSample(sampleId, sample);
    }

    public void addSampleFromSingleFile(String sampleId, double e, int l, int i, String filename, boolean nanook) {
        BlastSet sample = new BlastSet(e, l, i);
        sample.parseFile(filename, nanook);
        processNewSample(sampleId, sample);
    }
    
    public void compareGenes(boolean withMockRefs, String outname) {
        GeneSet allGenes = new GeneSet("All");
        GeneSet mockRefGenes = new GeneSet("MockRef");
        
        for (int i=0; i<geneSets.size(); i++) {
            allGenes.addGeneSet(geneSets.get(i));
        }
        
        if (withMockRefs) {
            cgp.parseListFile("/Users/leggettr/Desktop/BAMBI/mock_ref.txt", mockRefGenes);
            allGenes.addGeneSet(mockRefGenes);
        }

        try {
            PrintWriter pw = new PrintWriter(new FileWriter("/Users/leggettr/Desktop/BAMBI/" + outname));
       
            pw.print("Gene");
            
            for (int i=0; i<samples.size(); i++) {
                pw.printf("\t%s", sampleIds.get(i));
            }
            
            if (withMockRefs) {
                pw.printf("\tMockRef");
            }
            
            pw.print("\n");
            
            
            TreeMap<String, String> all = allGenes.getGenes();
            for (String gene : all.keySet()) {
                 pw.print(gene);                
                for (int i=0; i<geneSets.size(); i++) { 
                    GeneSet gs = geneSets.get(i);
                    int count = gs.getCount(gene);
                    //pw.printf("\t%s", count == 0 ? "NA":"1"); //Integer.toString(count));
                    pw.printf("\t%s", count == 0 ? "NA":Integer.toString(count));
                }
                
                if (withMockRefs) {
                    int count = mockRefGenes.getCount(gene);
                    //pw.printf("\t%s", count == 0 ? "NA":"1"); //Integer.toString(count));
                    pw.printf("\t%s", count == 0 ? "NA":Integer.toString(count));
                }
                
                pw.print("\n");
            }
            pw.close();
        } catch (Exception e) {
            System.out.println("Exception");
            e.printStackTrace();
        }                       
    }
    
    public void listAllGenes() {
        for (int i=0; i<geneSets.size(); i++) {
            String id = sampleIds.get(i);
            GeneSet gs = geneSets.get(i);
            System.out.println("SAMPLE "+id);
            gs.listGenes();
        }
    }
    
    public void writeSummaries() {
        for (int i=0; i<geneSets.size(); i++) {
            GeneSet gs = geneSets.get(i);
            gs.writeGeneCounts(cgp, "/Users/leggettr/Desktop/BAMBI/" + gs.getId() + "_all_genes.txt");
            //System.out.println("----> " + gs.getId());
            //gs.printGeneCounts(cgp);
        }
    }
    
    public static void processP8() {
        AMRAnalyser ma = new AMRAnalyser();
        ma.addSampleFromChunkSet("P8", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/BAMBI_1D_19092017/blastn_card", "all_Template_pass", "blastn_card");
        ma.addSampleFromSingleFile("IsolateMinION", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/blastn_nanopore_1.1.1.txt", false);        
        ma.addSampleFromSingleFile("IsolateIllumina", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/blastn_illumina_1.1.1.txt", false);
        ma.writeSummaries();
        ma.compareGenes(false, "test_P8_(200bp,80pc)_group(70pc).txt");        
    }
    
    public static void processP103M() {
        AMRAnalyser ma = new AMRAnalyser();
        ma.addSampleFromChunkSet("P103M_4pc", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P103M_P8Kpneu_BAMBI_10122018/blastn_card", "all_Template_pass", "blastn_card");
        ma.addSampleFromChunkSet("P103M_40pc", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P103M_P8Kpneu_BAMBI_19122018/blastn_card", "all_Template_pass", "blastn_card");
        ma.addSampleFromSingleFile("IsolateMinION", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/blastn_nanopore_1.1.1.txt", false);        
        ma.addSampleFromSingleFile("IsolateIllumina", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/blastn_nanopore_1.1.1.txt", false);
        ma.addSampleFromChunkSet("P103M_previous", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/20171220_1133_BAMBI_P103M_400ng_RAD4_20122017/blastn_card", "all_Template_pass", "blastn_card");
        ma.compareGenes(false, "test_P103M_all_(200bp,80pc)_group(70pc).txt");        
    }
    
    public static void processP10() {
        AMRAnalyser ma = new AMRAnalyser();
        ma.addSampleFromSingleFile("P10N_MinION", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/N79681_PRO977_P10n2_25082015_all_2D_pass.blastn.txt", true);   
        ma.addSampleFromSingleFile("P10N_Illumina", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/P10N_1m_1601_LIB18155_LDI15604_TAAGGCGA-AAGGAGTA_L001_R1_all.blastn.txt", true);           
        ma.addSampleFromSingleFile("P10R_MinION", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/N79681_PRO977_P10R1R3_2D_pass.blastn.txt", true);
        ma.addSampleFromSingleFile("P10R_Illumina", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/P10R_1m_1601_LIB18156_LDI15605_TAAGGCGA-CTAAGCCT_L001_R1_all.blastn.txt", true);
        ma.addSampleFromSingleFile("P10V_MinION", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/N79596_PRO977_P10V1V3_2D_pass.blastn.txt", true);
        ma.addSampleFromSingleFile("P10V_Illumina", 0.001, 200, 80, "/Users/leggettr/Desktop/BAMBI/P10/P10V_1601_LIB18157_LDI15606_CGTACTAG-GCGTAAGA_L001_R1_all.blastn.txt", true);
        ma.compareGenes(false, "test_P10_all_(200bp,80pc)_group(70pc).txt");        
        
        ma.listAllGenes();

        ma.writeSummaries();
    }
    
    public static void processFlongle() {
        AMRAnalyser ma = new AMRAnalyser();
        ma.addSampleFromChunkSet("MinION", 0.001, 200, 80, "/Users/leggettr/Desktop/Flongle/P129B_Flongle_Laptop/blastn_card",  "all_Template_pass", "blastn_card");
        ma.addSampleFromChunkSet("GridION", 0.001, 200, 80, "/Users/leggettr/Desktop/Flongle/P129B_Flongle_GridION/blastn_card", "all_Template_pass", "blastn_card");
        ma.compareGenes(false, "Flongle_(200bp,80pc)_group(70pc).txt");        
     }
    
    public static void main(String[] args) {
        CardGeneParser cgp = new CardGeneParser("/Users/leggettr/Desktop/BAMBI/all_genes_grouping_70pc.txt");
        BlastSet blastSampleCard = new BlastSet(0.001, 200, 80);
        BlastSet blastSampleCardB = new BlastSet(0.001, 200, 80);
        BlastSet blastIsolateMinionCard = new BlastSet(0.001, 200, 80);
        BlastSet blastIsolateIlluminaCard = new BlastSet(0.001, 200, 80);
        BlastSet blastSampleBacteria = null; 
        BlastSet blastPreviousCard = new BlastSet(0.001, 200, 80);
        BlastSet blastPreviousBacteria = null;
        GeneSet sampleAMRGenes = new GeneSet("Sample");
        GeneSet sampleBAMRGenes = new GeneSet("SampleB");
        GeneSet previousAMRGenes = new GeneSet("Previous");
        GeneSet isolateMinionAMRGenes = new GeneSet("IsolateMinION");
        GeneSet isolateIlluminaAMRGenes = new GeneSet("IsolateIllumina");
        GeneSet mockRefGenes = new GeneSet("MockRef");
        GeneSet allGenes = new GeneSet("All");
        boolean withMockRefs = false;
        boolean withPreviousSample = false;

        //processFlongle();
        processP10();
        System.out.println("Done!");
        System.exit(0);
        
        //blastSampleBacteria = new BlastSet();
        //blastSampleBacteria.parseReadToSpecies("/Users/leggettr/Desktop/BAMBI/all_reads_less_than_3000_subsample/blastn_bacteria", "all_Template_pass", "blastn_bacteria");
        
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/all_reads_less_than_3000_subsample/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);
        
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "MockWithSpike";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/all_reads_less_than_5000_subsample/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);

        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "MockWithSpike";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_subsample_lt3000_ss100000/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);

        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "MockWithSpikeT2";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_subsample_lt3000_ss100000_t2/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);
        
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "MockWithSpike50pc";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_50KbpMock_nbc_11032019_subsample_lt3000_ss100000/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);
        
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "MockInsilico";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //withMockRefs = true;
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_insilico_99999/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);

        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        String sampleName = "MockInsilico8";
        String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        withMockRefs = true;
        blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_insilico_88888/blastn_card", "all_Template_pass", "blastn_card");
        cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);        
        
        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //withPreviousSample = true;
        //String sampleName = "P103M_4pcSpike";
        String sampleNameB = "P103M_40pcSpike";

        //blastSampleBacteria.parseReadToSpecies("/Users/leggettr/Desktop/BAMBI/P103M_P8Kpneu_BAMBI_19122018/blastn_bacteria", "all_Template_pass", "blastn_bacteria");        

        //String sampleName = "P103M_4pc";
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/P103M_P8Kpneu_BAMBI_10122018/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);

        //String sampleName = "P103M_40pc";
        //blastSampleCardB.parseChunkSet("/Users/leggettr/Desktop/BAMBI/P103M_P8Kpneu_BAMBI_19122018/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCardB, sampleBAMRGenes, blastSampleBacteria);

        //blastPreviousCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/20171220_1133_BAMBI_P103M_400ng_RAD4_20122017/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastPreviousCard, previousAMRGenes, blastPreviousBacteria);
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";

        // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //String sampleName = "P8";
        //String outname = sampleName + "(200bp,80pc)_group(70pc).txt";
        //blastSampleCard.parseChunkSet("/Users/leggettr/Desktop/BAMBI/BAMBI_1D_19092017/blastn_card", "all_Template_pass", "blastn_card");
        //cgp.parseBlastAlignments(blastSampleCard, sampleAMRGenes, blastSampleBacteria);
        
        //blastIsolateMinionCard.parseFile("/Users/leggettr/Desktop/BAMBI/megablast_minion_card_1.1.1.txt", false);
        blastIsolateMinionCard.parseFile("/Users/leggettr/Desktop/BAMBI/blastn_nanopore_1.1.1.txt", false);
        cgp.parseBlastAlignments(blastIsolateMinionCard, isolateMinionAMRGenes, null);

        //blastIsolateIlluminaCard.parseFile("/Users/leggettr/Desktop/BAMBI/megablast_illumina_card_1.1.1.txt", false);
        blastIsolateIlluminaCard.parseFile("/Users/leggettr/Desktop/BAMBI/blastn_illumina_1.1.1.txt", false);
        cgp.parseBlastAlignments(blastIsolateIlluminaCard, isolateIlluminaAMRGenes, null);

        
        //System.out.println("\nMock...");
        //mockWithSpikeAMRGenes.listGenes();
        //System.out.println("\nIsolate MinION...");
        //isolateMinionAMRGenes.listGenes();
        //System.out.println("\nIsolate Illumina...");
        //isolateIlluminaAMRGenes.listGenes();
 
        allGenes.addGeneSet(sampleAMRGenes);
        allGenes.addGeneSet(isolateMinionAMRGenes);
        allGenes.addGeneSet(isolateIlluminaAMRGenes);
         
        if (withMockRefs) {
            cgp.parseListFile("/Users/leggettr/Desktop/BAMBI/mock_ref_only_good.txt", mockRefGenes);
            //System.out.println("\nMock ref:");
            //mockRefGenes.listGenes();
           allGenes.addGeneSet(mockRefGenes);
        }
        
        if (withPreviousSample) {
            allGenes.addGeneSet(previousAMRGenes);
            allGenes.addGeneSet(sampleBAMRGenes);
        }
                
        //System.out.println("\nALL:");
        //allGenes.listGenes();
               
        try {
            PrintWriter pw = new PrintWriter(new FileWriter("/Users/leggettr/Desktop/BAMBI/" + outname));
       
            if (withMockRefs) {
                //pw.printf("%70s %15s %15s %15s %11s %15s %s\n", "Gene", sampleName, "IsolateMinION", "IsolateIllumina", "MockRef", "Unexplained", "Bacteria");
                //pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", sampleName, "IsolateMinION", "IsolateIllumina", "MockRef", "Unexplained", "Bacteria");
                pw.printf("%s\t%s\t%s\t%s\t%s\n", "Gene", sampleName, "IsolateMinION", "IsolateIllumina", "MockRef");
            } else if (withPreviousSample) {
                pw.printf("%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", sampleName, sampleNameB, "IsolateMinION", "IsolateIllumina", "P103Previous");
            } else {                    
                pw.printf("%s\t%s\t%s\t%s\n", "Gene", sampleName, "IsolateMinION", "IsolateIllumina");
            }
            
            TreeMap<String, String> all = allGenes.getGenes();
            for (String gene : all.keySet()) {
                boolean sampleHasGene = sampleAMRGenes.hasGene(gene);
                boolean sampleBHasGene;
                boolean isolateMinIONHasGene = isolateMinionAMRGenes.hasGene(gene);
                boolean isolateIlluminaHasGene = isolateIlluminaAMRGenes.hasGene(gene);
                boolean mockRefHasGene = false;
                boolean previousSampleHasGene = false;
                boolean unexplained = sampleHasGene && !isolateMinIONHasGene && !isolateIlluminaHasGene && !mockRefHasGene;
                int sampleCount = sampleAMRGenes.getCount(gene);
                int sampleBCount = 0;
                int isolateMinIONCount = isolateMinionAMRGenes.getCount(gene);
                int isolateIlluminaCount = isolateIlluminaAMRGenes.getCount(gene);
                int mockCount = 0;
                int previousSampleCount = 0;
                boolean withCounts = false;
                
                if (withMockRefs) {
                    mockRefHasGene = mockRefGenes.hasGene(gene);
                    mockCount = mockRefGenes.getCount(gene);
                    unexplained = sampleHasGene && !isolateMinIONHasGene && !isolateIlluminaHasGene && !mockRefHasGene;
                } else if (withPreviousSample) {
                    sampleBHasGene = sampleBAMRGenes.hasGene(gene);
                    sampleBCount = sampleBAMRGenes.getCount(gene);
                    previousSampleHasGene = previousAMRGenes.hasGene(gene);
                    previousSampleCount = previousAMRGenes.getCount(gene);
                    unexplained = sampleHasGene && !isolateMinIONHasGene && !isolateIlluminaHasGene && !previousSampleHasGene;
                } else {
                    unexplained = sampleHasGene && !isolateMinIONHasGene && !isolateIlluminaHasGene;
                }
                
                if (withMockRefs) {
                    if (withCounts) {
                        pw.printf("%s\t%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":Integer.toString(sampleCount),
                                      isolateMinIONCount == 0 ? "NA":Integer.toString(isolateMinIONCount),
                                      isolateIlluminaCount == 0 ? "NA":Integer.toString(isolateIlluminaCount),
                                      mockCount == 0 ? "NA":Integer.toString(mockCount)
                                      //unexplained ? "1":"NA"
    //                                  sampleAMRGenes.getSpecies(gene)
                                      );
                    
                    } else {
                        //pw.printf("%70s %15s %15s %15s %15s %11s %s\n",
                        //              gene,
                        //              sampleHasGene ? "Y":" ",
                        //              isolateMinIONHasGene ? "Y":" ",
                        //              isolateIlluminaHasGene ? "Y":" ",
                        //              mockRefHasGene ? "Y":" ",
                        //              unexplained ? "*":"",
                        //              sampleAMRGenes.getSpecies(gene)
                        //              );
                        //pw.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        pw.printf("%s\t%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":"1",
                                      isolateMinIONCount == 0 ? "NA":"1",
                                      isolateIlluminaCount == 0 ? "NA":"1",
                                      mockCount == 0 ? "NA":"1"
    //                                  sampleAMRGenes.getSpecies(gene)
                                      );
                    }
                } else if (withPreviousSample) {
                    if (withCounts) {
                        pw.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":Integer.toString(sampleCount),
                                      sampleBCount == 0 ? "NA":Integer.toString(sampleBCount),
                                      isolateMinIONCount == 0 ? "NA":Integer.toString(isolateMinIONCount),
                                      isolateIlluminaCount == 0 ? "NA":Integer.toString(isolateIlluminaCount),
                                      previousSampleCount == 0 ? "NA":Integer.toString(previousSampleCount)
                                      );                        
                    } else {
                        pw.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":"1",
                                      sampleBCount == 0 ? "NA":"1",
                                      isolateMinIONCount == 0 ? "NA":"1",
                                      isolateIlluminaCount == 0 ? "NA":"1",
                                      previousSampleCount == 0 ? "NA":"1"
                                      );
                    }
                } else {
                    if (withCounts) {
                        pw.printf("%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":Integer.toString(sampleCount),
                                      isolateMinIONCount == 0 ? "NA":Integer.toString(isolateMinIONCount),
                                      isolateIlluminaCount == 0 ? "NA":Integer.toString(isolateIlluminaCount));                        
                    } else {
                        pw.printf("%s\t%s\t%s\t%s\n",
                                      gene + (unexplained ? "*":""),
                                      sampleCount == 0 ? "NA":"1",
                                      isolateMinIONCount == 0 ? "NA":"1",
                                      isolateIlluminaCount == 0 ? "NA":"1");                        
                        //pw.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
                        //              gene,
                        //              sampleHasGene ? "1":"0",
                        //              isolateMinIONHasGene ? "1":"0",
                        //              isolateIlluminaHasGene ? "1":"0",
                        //              unexplained ? "1":"0",
                        //              sampleAMRGenes.getSpecies(gene));
                    }
                }
            }
            pw.close();
        } catch (Exception e) {
            System.out.println("Exception");
            e.printStackTrace();
        }        
        
    }
    
}
