package org.abelhj.haplotect_utils;

import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.util.Collections;
import java.lang.Math;
import java.io.PrintStream;

import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.GenomeLoc;

public class HapCounter {
    
    private LinkedHashMap <String, ArrayList<BaseandQual[]> > haps=null;
    private LinkedHashMap <String, Double> popfreqs=null;
    private SnpPair pair=null;
    
    public HapCounter() {
        haps=new LinkedHashMap <String, ArrayList<BaseandQual[]> >();
    }
    
    public HapCounter(SnpPair pair) {

        this.pair=pair;
        haps=new LinkedHashMap <String, ArrayList<BaseandQual[]> >();
        for(String hh: pair.getFreqs().keySet()) {                          //only count haplotypes seen in reference popn
            haps.put(hh, new ArrayList<BaseandQual[]>());
        }
    }
    
    
    
    public void add(BaseandQual el1, BaseandQual el2) {
        char c1=(char) el1.getBase();
        char c2=(char) el2.getBase();
	BaseandQual[] els={el1, el2};
        String hh=c1+""+c2;
        if(haps.containsKey(hh)) {
            haps.get(hh).add(els);
        } 
    }

    public void print() {
      print(System.err);
    }
           
           
    public void print(PrintStream ps) {
        for(String ss: haps.keySet()) {
            print(ps, ss);
        }
    }
            


    public void print(PrintStream ps, String ss) {
        if(haps.containsKey(ss)) {
            for(BaseandQual[] plarray: haps.get(ss))
                ps.println(ss+"\t"+plarray[0]+"\t"+plarray[1]);
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
		Collections.reverse(counts);
        if(counts.size()>2)
			return counts.get(2);
        else
            return 0;
    }

    public int fourthHapCount() {                                            //number of reads matching fourth most common haplotype
        ArrayList<Integer> counts=new ArrayList<Integer>();
        for(String hh:haps.keySet()) 
            counts.add(haps.get(hh).size());
        Collections.sort(counts);
		Collections.reverse(counts);
        if(counts.size()>3)
			return counts.get(3);
        else
            return 0;
    }
    
    public void getPopFreqs() {
        popfreqs=new LinkedHashMap<String, Double>();
		int total=0;
		for(String hh: pair.getFreqs().keySet()) {
			total+=pair.getFreqs().get(hh); 
		}
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

    
    
    public double getLogLikMoreLogs(double alpha) {
        char[] nucs=new char[]{'A', 'C', 'G', 'T'};
        int[] errtype=new int[]{0,1,2,3};
        double result=0;
        getPopFreqs();
        LinkedHashMap<String, LinkedHashMap<String, LinkedHashMap<String, LinkedHashMap <String, Double> > > > lks=new LinkedHashMap<String, LinkedHashMap<String, LinkedHashMap<String, LinkedHashMap <String, Double> > > >();
        int hct=0;
        double maxlogprod=-999999;
        for(String h11:popfreqs.keySet()) {
            lks.put(h11, new  LinkedHashMap<String, LinkedHashMap<String, LinkedHashMap <String, Double> > >());
            for(String h12:popfreqs.keySet()) {
                lks.get(h11).put(h12, new LinkedHashMap<String, LinkedHashMap <String, Double> >());
                for(String h21:popfreqs.keySet()) {
                    lks.get(h11).get(h12).put(h21, new LinkedHashMap <String, Double>());
                    for(String h22:popfreqs.keySet()) {
                        
                        LinkedHashMap < Integer, LinkedHashMap<String, Double> >  probc=new LinkedHashMap < Integer, LinkedHashMap<String, Double> >();
                        LinkedHashMap < Integer, LinkedHashMap<String, Double> > probnc=new LinkedHashMap < Integer, LinkedHashMap<String, Double> >();
                        LinkedHashMap < Integer, Integer> countnc=new LinkedHashMap < Integer, Integer>();
                        LinkedHashMap < Integer, Integer> countc=new LinkedHashMap < Integer, Integer>();
                        
                        for(Integer err: errtype) {
                            
                            LinkedHashMap<String, Double> tempc=new LinkedHashMap<String, Double>();
                            LinkedHashMap<String, Double> tempnc=new LinkedHashMap<String, Double>();
                            
                            int tempcountc=0;
                            int tempcountnc=0;
                            
                            for(Character c: nucs) {
                                for(Character d : nucs) {
                                    tempc.put(c+""+d, 0.0);
                                    tempnc.put(c+""+d, 0.0);
                                }
                            }
                            
                            for(String ss: tempc.keySet()) {
                                if(scoreMatch(ss, h11)==err) {
                                    tempcountnc++;
                                    tempnc.put(ss, tempnc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h12)==err) {
                                    tempcountnc++;
                                    tempnc.put(ss, tempnc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h21)==err) {
                                    tempcountc++;
                                    tempc.put(ss, tempc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h22)==err) {
                                    tempcountc++;
                                    tempc.put(ss, tempc.get(ss)+1);
                                }
                            }
                            probc.put(err, tempc);
                            probnc.put(err, tempnc);
                            countc.put(err, tempcountc);
                            countnc.put(err, tempcountnc);
                            
                        }
                        
                        double logprod=0;
                        for(String hh: haps.keySet()) {
                            ArrayList<BaseandQual[]> els=haps.get(hh);
                            for(BaseandQual[] pel: els) {
                                double sum=0.0;
                                for(Integer err: errtype) {
                                    
                                    sum+=(1-alpha)*probnc.get(err).get(hh)/countnc.get(err)+alpha*probc.get(err).get(hh)/countc.get(err);
                                    sum*=errProb(err, pel);
                                }
                                logprod+=Math.log(sum);
                            }
                        }
                        if(hct==0) {
                            maxlogprod=logprod;
                        } else if(logprod>maxlogprod) {
                            maxlogprod=logprod;
                        }
                        hct++;
                        lks.get(h11).get(h12).get(h21).put(h22, logprod);
                    }
                }
            }
        }
        double logsum=0;
        for(String h11:popfreqs.keySet()) {
            for(String h12: popfreqs.keySet()) {
                for(String h21:popfreqs.keySet()) {
                    for(String h22: popfreqs.keySet()) {
                        logsum+=Math.exp(lks.get(h11).get(h12).get(h21).get(h22)-maxlogprod);
                    }
                }
            }
        }
        logsum=Math.log(logsum)+maxlogprod;
        return logsum;
    }
    
    
    public double getLogLik(double alpha) {
        char[] nucs={'A', 'C', 'G', 'T'};
        int[] errtype={0,1,2,3};
        double result=0;
        getPopFreqs();
        for(String h11:popfreqs.keySet()) {
            for(String h12:popfreqs.keySet()) {
                for(String h21:popfreqs.keySet()) {
                    for(String h22:popfreqs.keySet()) {
                        
                        LinkedHashMap < Integer, LinkedHashMap<String, Double> >  probc=new LinkedHashMap < Integer, LinkedHashMap<String, Double> >();
                        LinkedHashMap < Integer, LinkedHashMap<String, Double> > probnc=new LinkedHashMap < Integer, LinkedHashMap<String, Double> >();
                        LinkedHashMap < Integer, Integer> countnc=new LinkedHashMap < Integer, Integer>();
                        LinkedHashMap < Integer, Integer> countc=new LinkedHashMap < Integer, Integer>();

                        
                        
                        for(Integer err: errtype) {
                            
                            LinkedHashMap<String, Double> tempc=new LinkedHashMap<String, Double>();
                            LinkedHashMap<String, Double> tempnc=new LinkedHashMap<String, Double>();
                            
                            int tempcountc=0;
                            int tempcountnc=0;
                            
                            for(Character c: nucs) {
                                for(Character d : nucs) {
                                    tempc.put(c+""+d, 0.0);
                                    tempnc.put(c+""+d, 0.0);
                                }
                            }
                           
                            for(String ss: tempc.keySet()) {
                                if(scoreMatch(ss, h11)==err) {
                                    tempcountnc++;
                                    tempnc.put(ss, tempnc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h12)==err) {
                                    tempcountnc++;
                                    tempnc.put(ss, tempnc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h21)==err) {
                                    tempcountc++;
                                    tempc.put(ss, tempc.get(ss)+1);
                                }
                                if(scoreMatch(ss, h22)==err) {
                                    tempcountc++;
                                    tempc.put(ss, tempc.get(ss)+1);
                                }
                            }
                            probc.put(err, tempc);
                            probnc.put(err, tempnc);
                            countc.put(err, tempcountc);
                            countnc.put(err, tempcountnc);
                            
                        }
                        
                        double prod=1.0;
                        for(String hh: haps.keySet()) {
                            ArrayList<BaseandQual[]> els=haps.get(hh);
                            for(BaseandQual[] pel: els) {
                                double sum=0.0;
                                for(Integer err: errtype) {
                                    sum+=(1-alpha)*probnc.get(err).get(hh)/countnc.get(err)+alpha*probc.get(err).get(hh)/countc.get(err);
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
        return Math.log(result);
    }


    public static int scoreMatch(String hh, String kk) {
        int score=0;
        score+=(hh.charAt(0)==kk.charAt(0)?1:0);
        score+=2*(hh.charAt(1)==kk.charAt(1)?1:0);
        return score;
    }
    
    public static double errProb(int errtype, BaseandQual[] el) {
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
            ret+=hh+":"+haps.get(hh).size()+";";
        }
        return ret;
    }
    
                          
         
}