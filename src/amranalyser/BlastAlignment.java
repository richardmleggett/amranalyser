package amranalyser;

import java.util.ArrayList;

public class BlastAlignment {
    // qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
    private String qseqid;
    private String sseqid;
    private float pident;
    private int length;
    private int mismatch;
    private int gapopen;
    private int qstart;
    private int qend;
    private int sstart;
    private int send;
    private double evalue;
    private double bitscore;
    private String stitle;
    private int taxonId = -1;
    private boolean validAlignment = false;
    private long leafNode = 0;
    private ArrayList<Long> taxonIdPath;
    
    public BlastAlignment(String a, boolean nanook, double maxEValue, int minLength, int minIdent) {
        String fields[] = a.split("\t", 14);
        boolean debug = false;
        
        // non-nanook: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle
        
        if (fields.length >= 13) {
            try {
                qseqid = fields[0];
                
                if (qseqid.equals("D3NG5HQ1:381:HCFH2BCXX:1:2204:7128:37328")) {
                    debug = true;
                }                
                
                sseqid = fields[1];
                pident = Float.parseFloat(fields[2]);
                length = Integer.parseInt(fields[3]);
                mismatch = Integer.parseInt(fields[4]);
                gapopen = Integer.parseInt(fields[5]);
                qstart = Integer.parseInt(fields[6]);
                qend = Integer.parseInt(fields[7]);
                sstart = Integer.parseInt(fields[8]);
                send = Integer.parseInt(fields[9]);
                evalue = Double.parseDouble(fields[10]);
                bitscore = Double.parseDouble(fields[11]);
                
                if (nanook) {                
                    stitle = fields[12];
                    validAlignment = true;
                    if (fields.length == 14) {
                        if (!fields[13].equals("N/A")) {
                            taxonId = Integer.parseInt(fields[13]);
                        }
                    } else {
                        //System.out.println("No taxonid");
                    }
                } else {
                    //slen = Integer.parseInt(fields[12]);
                    stitle = fields[13];
                    validAlignment = true;
                }                
            } catch (Exception e) {
                System.out.println("Error parsing alignment - incomplete file?");
                e.printStackTrace();
                System.exit(1);
            }
        } else {
            System.out.println("Unknown alignment - incomplete file? ["+a+"]");
        }
        
        if (validAlignment) {
            if (evalue > maxEValue) {
                if (debug) System.out.println("Alignment rejected for evalue ("+evalue+") - "+a);
                validAlignment = false;
            }
            
            // 50 200
            if (length < minLength) {
                if (debug) System.out.println("Alignment rejected for length ("+length+") - "+a);
                validAlignment = false;
            }
            
            // 75 60 80
            if (pident < minIdent) {
                if (debug) System.out.println("Alignment rejected for pident ("+pident+") - "+a);
                validAlignment = false;
            }
        }
        
        if ((debug) && (validAlignment)) {
            System.out.println("Accepted: "+a);
        }
                
    }
    
    public boolean isValidAlignment() {
        return validAlignment;
    }
    
    public String getQueryId() {
        return qseqid;
    }
    
    public double getEValue() {
        return evalue;
    }
    
    public float getPercentIdentity() {
        return pident;
    }
    
    public String getSubjectId() {
        return sseqid;
    }
    
    public String getSubjectTitle() {
        return stitle;
    }
    
    public Double getBitScore() {
        return bitscore;
    }
    
    public int getLength() {
        return length;
    }
    
//    public void cacheTaxonIdPath(Taxonomy t) {
//        if (taxonId == -1) {
//            Long id = t.parseTaxonomyToId(stitle);
//
//            if (id != null) {
//                leafNode = id;
//            } else {
//                leafNode = 0;
//            }
//        } else {
//            leafNode = taxonId;
//        }
//
//        taxonIdPath = t.getTaxonIdPathFromId(leafNode);
//    }
//    
//    public ArrayList getTaxonIdPath() {
//        return taxonIdPath;
//    }
//    
//    public int getTaxonLevel() {
//        if (taxonIdPath != null) {
//            return taxonIdPath.size(); // 1-offset
//        }
//        return 0;
//    }
//    
//    public long getLeafNode() {
//        return leafNode;
//    }
//    
//    // Note level is 1-offset
//    public long getTaxonNode(int level) {
//        if (level <= taxonIdPath.size()) {
//            return taxonIdPath.get(taxonIdPath.size() - level);
//        }
//        return 0;
//    }
    
    public int getQueryStart() {
        return qstart;
    }

    public int getQueryEnd() {
        return qend;
    }

    public int getSubjectStart() {
        return sstart;
    }

    public int getSubjectEnd() {
        return send;
    }

}
