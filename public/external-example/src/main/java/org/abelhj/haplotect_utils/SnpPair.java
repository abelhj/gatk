package org.abelhj.haplotect_utils;

import org.broadinstitute.gatk.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.LinkedHashMap;

public class SnpPair {

    private GenomeLoc snp1=null;
    private GenomeLoc snp2=null;
    private String info=null;
    private LinkedHashMap<String, Integer> hapFreq=null;

    public SnpPair() {}

    public SnpPair(GenomeLoc v1, GenomeLoc v2) {
        snp1=v1;
        snp2=v2;
    }

    public SnpPair(GenomeLoc v1, GenomeLoc v2, String ss) {
        this(v1, v2, ss, false);    
    }

    public SnpPair(GenomeLoc v1, GenomeLoc v2, String ss, boolean unif) {
	snp1=v1;
	snp2=v2;
	info=ss;
	hapFreq=new LinkedHashMap<String, Integer>();
	String[] spl=info.split("\\t");
	if(!unif) {
	    String[] freqs=spl[7].split(";");
	    for(int i=0; i<freqs.length; i++) {
		String[] kv=freqs[i].split(":");
		hapFreq.put(kv[0], Integer.parseInt(kv[1]));
	    }
	} else {
	    String str=(new StringBuilder()).append(spl[3]).append(spl[5]).toString();
	    hapFreq.put(str, 100);
	    str=(new StringBuilder()).append(spl[3]).append(spl[6]).toString();
	    hapFreq.put(str, 100);
	    str=(new StringBuilder()).append(spl[4]).append(spl[5]).toString();
            hapFreq.put(str, 100);
	    str=(new StringBuilder()).append(spl[4]).append(spl[6]).toString();
            hapFreq.put(str, 100);
	}
    }
    
    public LinkedHashMap<String, Integer> getFreqs() {
        return hapFreq;
    }
    
    public int distance() {
        return snp2.distance(snp1);
    }
    
    public String pairInfo() {
        return info;
    }

    public ArrayList<GenomeLoc> getSnps() {
        ArrayList<GenomeLoc> snps=new ArrayList<GenomeLoc>();
        snps.add(snp1);
        snps.add(snp2);
        return snps;
    }
    
    public boolean matches(GenomeLoc loc) {
        if(loc.equals(snp1) || loc.equals(snp2))
            return true;
        else return false;
    }
    
    public int whichSnp(GenomeLoc loc) {
        if(loc.equals(snp1))
            return 1;
        else if ( loc.equals(snp2))
            return 2;
        else return -1;
    }


    public String toString() {
        String ret="";
        ret+="["+snp1+","+snp2+"]";
        return ret;
    }
}
	