/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package amranalyser;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.TreeMap;

/**
 *
 * @author leggettr
 */
public class GeneSet {
    String sampleId;
    private TreeMap<String, String> genesToSpecies = new TreeMap<String, String>(); 
    private TreeMap<String, Integer> genesToCount = new TreeMap<String, Integer>();
    private TreeMap<String, Integer> aroCounts = new TreeMap<String, Integer>();
    
//new TreeMap<>(new Comparator<String>() {
//            @Override
//            public int compare(String s1, String s2) {
//                return s2.compareTo(s1);
//            }
//        });
    
    public GeneSet(String id) {
        sampleId = id;
    }
    
    public String getId() {
        return sampleId;
    }
    
    public void addGeneSet(GeneSet g) {
        TreeMap<String, String> genesToAdd = g.getGenes();
        TreeMap<String, Integer> geneCountsToAdd = g.getGeneCounts();
        
        for (String gene : genesToAdd.keySet()) {
            String c = genesToAdd.get(gene);
            int count = geneCountsToAdd.get(gene);
            if (genesToSpecies.containsKey(gene)) {
                c = c + "," + genesToSpecies.get(gene);
                count = count + genesToCount.get(gene);
            }
            
            genesToSpecies.put(gene, c);
            genesToCount.put(gene, count);
        }        
    }
    
    public void addGene(String a, String g, String gg, String b) {
        if (aroCounts.containsKey(a)) {
            int count = aroCounts.get(a);
            count +=1;
            aroCounts.put(a, count);
        } else {
            aroCounts.put(a, 1);
        }
        
        if (genesToSpecies.containsKey(gg)) {
            String c = genesToSpecies.get(gg);
            int count = genesToCount.get(gg);
            c = c + "," + b;
            count += 1;
            genesToSpecies.put(gg, c);
            genesToCount.put(gg, count);
        } else {
            genesToSpecies.put(gg, b);
            genesToCount.put(gg, 1);
        }
    }
    
    public void printGeneCounts(CardGeneParser cgp) {
        for (String aro : aroCounts.keySet()) {
            System.out.println(aroCounts.get(aro) + "\t" + aro + "\t" + cgp.geneToGroup(aro));
        }
    }
    
    public void writeGeneCounts(CardGeneParser cgp, String filename) {
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(filename));
            for (String aro : aroCounts.keySet()) {
                String gene = cgp.aroToGene(aro);
                String geneGroup = cgp.geneToGroup(gene);
                pw.println(aroCounts.get(aro) + "\t" + aro + "\t" + gene + "\t" + geneGroup);
            }
            pw.close();
        } catch (Exception e) {
            System.out.println("Exception");
            e.printStackTrace();
        }  
    }
    
    public void listGenes() {
        for (String gene : genesToSpecies.keySet()) {
            System.out.println(gene);
            //System.out.println(gene + "(" + genes.get(gene) + ")");
        }
    }
    
    public TreeMap<String, String> getGenes() {
        return genesToSpecies;
    }
    
    public TreeMap<String, Integer> getGeneCounts() {
        return genesToCount;
    }
    
    public boolean hasGene(String g) {
        if (genesToSpecies.containsKey(g)) {
            return true;
        } else {
            return false;
        }
    }
    
    public String getSpecies(String g) {
        if (genesToSpecies.containsKey(g)) {
            return genesToSpecies.get(g);
        } else {
            return "";
        }
    }

    public int getCount(String g) {
        if (genesToCount.containsKey(g)) {
            return genesToCount.get(g);
        } else {
            return 0;
        }
    }
}
