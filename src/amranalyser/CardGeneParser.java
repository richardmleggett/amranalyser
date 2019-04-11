/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package amranalyser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

/**
 *
 * @author leggettr
 */
public class CardGeneParser {
    private AROMap aroMap = new AROMap();
    private HashMap<String, String> geneToGroup = new HashMap<String, String>(); 

    
    public CardGeneParser(String geneGroupsFile) {
        AROMap.readMapFile("/Users/leggettr/Documents/Databases/CARD_1.1.1_Download_17Oct16/aro.csv");
        readGeneGroups(geneGroupsFile);
    }
    
    public void readGeneGroups(String filename) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                String groupName = "group " + fields[0];
                int nGenes = Integer.parseInt(fields[1]);
                if (nGenes > 1) {
                    for (int i=0; i<nGenes; i++) {
                        if (geneToGroup.containsKey(fields[2+i])) {
                            if (!geneToGroup.get(fields[2+i]).equals(groupName)) {
                                System.out.println("Error: " + fields[2+i] + " alreaday linked to " + geneToGroup.get(fields[2+i]) + " - can't link to " + groupName);
                                System.exit(1);
                            }
                        } else {
                            geneToGroup.put(fields[2+i], groupName);
                        }
                    }
                }
            }
                    
            br.close();
        } catch (Exception e) {
            System.out.println("Exception");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public String aroToGene(String aro) {
        String gene = aroMap.getNameFromAccession(aro);  
        return gene;
    }
    
    public String geneToGroup(String gene) {
        if (geneToGroup.containsKey(gene)) {
            return geneToGroup.get(gene);
        } else {
            return gene;
        }
    }
    
    public void parseBlastAlignments(BlastSet bs, GeneSet gs, BlastSet bsLookup) {
        for (int i=0; i<bs.getNumberOfAlignments(); i++) {
            BlastAlignment ba = bs.getAlignment(i);
            String subjId = ba.getSubjectId();
            int aroPos = subjId.indexOf("|ARO:");
            String aro = "";

            if (aroPos > 0) { 
                String aroStart = subjId.substring(aroPos+1);
                aro = aroStart.substring(0, aroStart.indexOf('|'));
            } else {
                System.out.println("WARNING: Couldn't find ARO in " + subjId);
                aro = subjId;
            }
            
            String gene = aroToGene(aro);
            String geneGroup = geneToGroup(gene);
            String bacteria = ba.getQueryId();
            //System.out.println(ba.getQueryId() + " maps to " + gene);
            
            if (bsLookup == null) {
                bacteria = ba.getQueryId();
            } else {
                //System.out.println("Looking up "+ba.getQueryId());
                bacteria = "\"" + bsLookup.getSpeciesForRead(ba.getQueryId()) + "\"";
            }
            
            gs.addGene(aro, gene, geneGroup, bacteria);
        }
    }
    
    public void parseListFile(String filename, GeneSet gs) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = br.readLine()) != null) {

                if (line.length() > 5) {
                    int aroPos = line.indexOf("|ARO:");
                    String aro = "";

                    if (aroPos > 0) { 
                        String aroStart = line.substring(aroPos+1);
                        aro = aroStart.substring(0, aroStart.indexOf('|'));
                    } else {
                        System.out.println("WARNING: Couldn't find ARO in " + line);
                        aro = line;   
                    }
                    
                    String gene = aroToGene(aro);
                    String geneGroup = geneToGroup(gene);
                    gs.addGene(aro, gene, geneGroup, "hit");
                }
            }
                    
            br.close();
        } catch (Exception e) {
            System.out.println("Exception");
            e.printStackTrace();
            System.exit(1);
        }
        
    }
}
