package org.abelhj;

import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.util.Collections;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

public class HapCounter {
    
    private LinkedHashMap <String, ArrayList<PileupElement[]> > haps=null;
    private LinkedHashMap <String, Double> popfreqs=null;
    private SnpPair pair=null;
    
    public HapCounter() {
        haps=new LinkedHashMap <String, ArrayList<PileupElement[]> >();
    }
    
    public HapCounter(SnpPair pair) {

        this.pair=pair;
        haps=new LinkedHashMap <String, ArrayList<PileupElement[]> >();
        for(String hh: pair.getFreqs().keySet()) {                          //only count haplotypes seen in reference popn
            haps.put(hh, new ArrayList<PileupElement[]>());
        }
    }
    
    
    
    public void add(PileupElement el1, PileupElement el2) {
        char c1=(char) el1.getBase();
        char c2=(char) el2.getBase();
        PileupElement[] els=new PileupElement[2];
        els[0]=el1;
        els[1]=el2;
        String hh=c1+""+c2;
        if(haps.containsKey(hh)) {
            haps.get(hh).add(els);
        }
    }
           
           
    public void print() {
        for(String ss: haps.keySet()) {
            print(ss);
        }
    }
            


    public void print(String ss) {
        if(haps.containsKey(ss)) {
            for(PileupElement[] plarray: haps.get(ss))
                System.out.println(ss+"\t"+plarray[0]+"\t"+plarray[1]);
        }
    }
    
    public int totalCount() {                                               //total number of reads spanning snp pair (that match haplotype in ref popn)
        int count=0;
        for(String ss: haps.keySet())
            count+=haps.get(ss).size();
        return count;
    }


    
    public int numUniqueObsHaps() {
        int count=0;
        for(String ss: haps.keySet()) {
            if(haps.get(ss).size()>0)
                count++;
        }
        return count;
    }
           
    
    public int thirdHapCount() {                                            //number of reads matching third most common haplotype
        ArrayList<Integer> counts=new ArrayList<Integer>();
        for(String hh:haps.keySet()) 
            counts.add(haps.get(hh).size());
        Collections.sort(counts);
        if(counts.size()>2)
            return counts.get(counts.size()-3);
        else
            return 0;
    }
    
    public void getPopFreqs() {
        popfreqs=new LinkedHashMap<String, Double>();
        int total=totalCount();
        for(String hh: pair.getFreqs().keySet()) {
            int ct=pair.getFreqs().get(hh);
            double ff=-1;
            if(ct==0) 
                ff=0.005;
            else if (ct==total) 
                ff=0.995;
            else 
                ff=1.0*ct/total;
            popfreqs.put(hh, ff);
        }
    }
    
    public double getLogLik(double alpha) {
        char[] nucs=new char[]{'A', 'C', 'G', 'T'};
        int[] errtype=new int[]{0,1,2,3};
        double result=0;
        getPopFreqs();
        for(String h11:popfreqs.keySet()) {
            for(String h12:popfreqs.keySet()) {
                for(String h21:popfreqs.keySet()) {
                    for(String h22:popfreqs.keySet()) {
                        double prod=1.0;
                        for(String hh: haps.keySet()) {
                            ArrayList<PileupElement[]> els=haps.get(hh);
                            for(PileupElement[] pel: els) {
                                double sum=0.0;
                                for(Integer err: errtype) {
                                    LinkedHashMap<String, Double> probc=new LinkedHashMap<String, Double>();
                                    LinkedHashMap<String, Double> probnc=new LinkedHashMap<String, Double>();
                                    for(Character c: nucs) {
                                        for(Character d : nucs) {
                                            probc.put(c+""+d, 0.0);
                                            probnc.put(c+""+d, 0.0);
                                        }
                                    }
                                    int countc=0;
                                    int countnc=0;
                                    for(String ss: probc.keySet()) {
                                        if(scoreMatch(ss, h11)==err) {
                                            countnc++;
                                            probnc.put(ss, probnc.get(ss)+1);
                                        }
                                        if(scoreMatch(ss, h12)==err) {
                                            countnc++;
                                            probnc.put(ss, probnc.get(ss)+1);
                                        }
                                        if(scoreMatch(ss, h21)==err) {
                                            countc++;
                                            probc.put(ss, probc.get(ss)+1);
                                        }
                                        if(scoreMatch(ss, h22)==err) {
                                            countc++;
                                            probc.put(ss, probc.get(ss)+1);
                                        }
                                    }
				       
                                    sum+=(1-alpha)*probnc.get(hh)/countnc+alpha*probc.get(hh)/countc;
                                    sum*=errProb(err, pel);
                                }
                                prod*=sum;
                            }
                        }
                        result+=prod;
                    }
                }
            }
        }
        return result;
    }

    public static int scoreMatch(String hh, String kk) {
        int score=0;
        score+=(hh.charAt(0)==kk.charAt(0)?1:0);
        score+=2*(hh.charAt(1)==kk.charAt(1)?1:0);
        return score;
    }
    
    public static double errProb(int errtype, PileupElement[] el) {
        double q1=Math.pow(10, el[0].getQual()/-10.0);
        double q2=Math.pow(10, el[1].getQual()/-10.0);
        double ret=-1;
        if(errtype==0) 
            ret=q1*q2;
        else if (errtype==1) 
            ret=(1-q1)*q2;
        else if (errtype==2)
            ret=q1*(1-q2);
        else if (errtype==3)
            ret=(1-q1)*(1-q2);
        else {
            System.err.println("bad error type"+errtype);
            System.exit(1);
        }
        return ret;
    }
    
        
    
    public String countString() {
        String ret="";
        for(String hh : haps.keySet()) {
            ret+=hh+"\t"+haps.get(hh).size()+"\t";
        }
        return ret;
    }
    
                          
         
}