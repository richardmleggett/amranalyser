/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package amranalyser;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author leggettr
 */
public class BlastSet {
    private ArrayList<BlastAlignment> alignments = new ArrayList<BlastAlignment>();
    private HashMap<String, String> readIdToSpecies = new HashMap<String, String>(); 
    private double maxEValue = 0.001;
    private int minLength = 2000000;
    private int minIdent = 200;
            
    public BlastSet(double e, int l, int i) {
        maxEValue = e;
        minLength = l;
        minIdent = i;
    }
    
    public void parseChunkSet(String directory, String prefix, String midfix) {
        int c = 0;
        boolean foundNext = true;
        
        while((foundNext) && (c < 200)) {
            String filename = directory + File.separator + prefix + "_" + c + "_" + midfix + ".txt";
            File f = new File(filename);
            if (f.exists()) {
                System.out.println("Loading " + filename);
                BlastFile bf = new BlastFile(filename, true, maxEValue, minLength, minIdent);
                bf.parseAlignments(alignments);
                c++;
            } else {
                foundNext = false;
            }
        }            
    }
    
    public void parseReadToSpecies(String directory, String prefix, String midfix) {
        int c = 0;
        boolean foundNext = true;
        
        while(foundNext) {
            String filename = directory + File.separator + prefix + "_" + c + "_" + midfix + ".txt";
            File f = new File(filename);
            if (f.exists()) {
                System.out.println("Loading " + filename);
                BlastFile bf = new BlastFile(filename, true, maxEValue, minLength, minIdent);
                bf.parseIdToSpecies(readIdToSpecies);
                c++;
            } else {
                foundNext = false;
            }
        }            
    }
    
    public String getSpeciesForRead(String r) {
        if (readIdToSpecies.containsKey(r)) {
            return readIdToSpecies.get(r);
        } else {
            return r;            
        }
    }
    
    public void parseFile(String filename, boolean nanook) {
        System.out.println("Loading " + filename);
        BlastFile bf = new BlastFile(filename, nanook, maxEValue, minLength, minIdent);
        bf.parseAlignments(alignments);        
    }
    
    public int getNumberOfAlignments() {
        return alignments.size();
    }
    
    public BlastAlignment getAlignment(int i) {
        if (i < alignments.size()) {
            return alignments.get(i);
        } else {
            return null;
        }
    }
}
