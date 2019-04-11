/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package amranalyser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author leggettr
 */
public class BlastFile {
    String filename;
    boolean nanookGenerated = true;    
    private double maxEValue = 0.001;
    private int minLength = 2000000;
    private int minIdent = 200;
    
    public BlastFile(String f, boolean n, double e, int l, int i) {
        filename = f;        
        nanookGenerated = n;
        maxEValue = e;
        minLength = l;
        minIdent = i;
    }
    
    public void parseAlignments(ArrayList<BlastAlignment> alignments) {
        String lastId = "";
        //Double lastE = 0.0;
        CoordinateList cl = new CoordinateList();
        int lineNo = 0;
        boolean debug = false;

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.length() > 0) {
                    BlastAlignment ba = new BlastAlignment(line, nanookGenerated, maxEValue, minLength, minIdent);
                    lineNo++;
                    
                    if (ba.getQueryId().equals("30e029d6-260b-4b79-ad5d-0bc243b6573d_Basecall_2D_2d")) {
                        System.out.println("DEBUG");
                        debug = true;
                    } else {
                        debug = false;
                    }
                    
                    if (ba.isValidAlignment()) {
                        if (debug) System.out.println("Valid alignment");
                        if (ba.getQueryId().equals(lastId)) {
                            int overlap = cl.getOverlap(ba.getQueryStart(), ba.getQueryEnd());
                            if (overlap < (ba.getLength() / 10)) {
                                //System.out.println("Accepting overlap at line "+lineNo+" of size "+overlap);
                                cl.add(ba.getQueryStart(), ba.getQueryEnd());
                                alignments.add(ba);                        
                            } else {
                                //System.out.println("Overlaps at line "+lineNo+ " of size "+overlap);
                            }
                        } else {
                            //System.out.println(lastId + " count " + cl.getCount());
                            cl = new CoordinateList();
                            alignments.add(ba);
                            cl.add(ba.getQueryStart(), ba.getQueryEnd());
                        }

                        //cl.add(ba.getQueryStart(), ba.getQueryEnd());
                        lastId = ba.getQueryId();
                        //lastE = ba.getEValue();
                    } else {
                        if (debug) System.out.println("Not Valid alignment: "+line);
                    }
                } 
            }
                    
            br.close();
        } catch (Exception e) {
            System.out.println("parseAlignment exception");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    public void parseIdToSpecies(HashMap<String, String> species) {
        String lastId = "";

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = br.readLine()) != null) {
                BlastAlignment ba = new BlastAlignment(line, nanookGenerated, maxEValue, minLength, minIdent);

                if (ba.isValidAlignment()) {                    
                    // Take top hit for now
                    if (!ba.getQueryId().equals(lastId)) {
                        species.put(ba.getQueryId(), ba.getSubjectTitle());
                    }
                    
                    lastId = ba.getQueryId();
                }
            }
                    
            br.close();
        } catch (Exception e) {
            System.out.println("parseAlignment exception");
            e.printStackTrace();
            System.exit(1);
        }

    }
    
}
