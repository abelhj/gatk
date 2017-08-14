package org.abelhj.utils;

import java.util.List;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.PrintStream;

import org.broadinstitute.gatk.utils.GenomeLoc;


import htsjdk.samtools.SAMFlag;

public class BaseFlagMapAmpBC {

    private Map<String, Map<String, BaseFlagMap> > map=null;
    private Map<String, Map<String, BaseFlagMap> > bcmap=null;
    private int minBCPerAmp=50;
    private char refBase='\0';
    private double minPercentRG=0.8;
    private int minCountPerBC=1;
    private GenomeLoc loc=null;
    private double diffVAF=-99;
    private double pvalUncorr=-99;
    private int refPos;
    private int namplicons=-99;
    private String VAFString=null;
    private double maxVAF=-99;

    public BaseFlagMapBC(GenomeLoc pos, char ref, int minPerAmp, double minPct, int minCount){
	refBase=ref;
	map=new LinkedHashMap<Integer, Map<String, BaseFlagMap> >();
	bcmap=new LinkedHashMap<String, Map<Integer, BaseFlagMap> >();
	endPos=new LinkedHashMap<Integer, Integer>();
	loc=pos;
	refPos=loc.getStart();
	minBCPerAmp=minPerAmp;
	minPercentRG=minPct;
	minCountPerBC=minCount;
    }
    
    public void add(BaseFlagAmpBC bfl) {

	String amp=bfl.getAmp();
	String barcode=bfl.getBarcode();

	if(map.containsKey(amp) && map.get(amp).containsKey(barcode)) {
	    map.get(amp).get(barcode).add(bfl);
            bcmap.get(barcode).get(amp).add(bfl);
	} else {
	    BaseFlagMap bfm=new BaseFlagMap();
	    bfm.add(bfl);

	    if(!map.containsKey(amp)) {
		map.put(amp, new LinkedHashMap<String, BaseFlagMap>());
	    } 
	    if(!map.get(amp).containsKey(barcode)) {
		map.get(amp).put(barcode, bfm);
	    }

	    if(!bcmap.containsKey(barcode)) {
		bcmap.put(barcode, new LinkedHashMap<String, BaseFlagMap>());
	    }
	    if(!bcmap.get(barcode).containsKey(amp)) {
		bcmap.get(barcode).put(amp, bfm);
	    }
	}
    }


    public void calcAmpBias(boolean debug, char altBase) {

	List<Integer> refct=new ArrayList<Integer>();
	List<Integer> altct=new ArrayList<Integer>();
	List<String> amps=new ArrayList<String>();
	String str="";

	for (String amp : map.keySet()) {

	    if(map.get(amp).keySet().size()>=minBCPerAmp) {
       
		BaseFlagMap bf=new BaseFlagMap();
		for(String bc : map.get(amp).keySet()) {
		    
		    BaseFlagMap bfm=map.get(amp).get(bc);
		    char maxBase=map.get(amp).get(bc).maxBase(false);
		    int maxCt=bfm.sum(maxBase);
		    int total=bfm.sum();
		    if( maxCt*1.0/total > minPercentRG && total >= minCountPerBC) {
			bf.add(new BaseFlag(maxBase, SAMFlag.FIRST_OF_PAIR.intValue()));
		    }		
		}
		int altct_temp=bf.sum(altBase);
		int refct_temp=bf.sum(refBase);
		refct.add(refct_temp);
		altct.add(altct_temp);
		amps.add(amp);
		str+=(amp+":"+refct_temp+","+altct_temp+";");
	    }
	}
	diffVAF=ChiSquareUtils.maxDiffVAF(refct, altct);
	pvalUncorr=ChiSquareUtils.chiSquare(refct, altct);
	maxVAF=ChiSquareUtils.maxVAF(refct, altct);
        VAFString=str;
	namplicons=refct.size();
	if(debug) {
	  System.err.println(diffVAF+"\t"+pvalUncorr);
	  for(int i=0; i<refct.size(); i++) {
	      System.err.println(i+"\t"+refct.get(i)+"\t"+altct.get(i)+"\t"+amps.get(i));
	  }
	}
    }


    public void valueUnsetError(String var) {
	System.err.println("Error:  Value of "+var+"unset--need to calculate amplicon bias first.\n");
	System.exit(1);
    }

    public String getVAFString() {
	if(VAFString==null) {
	    valueUnsetError("VAF string");
	}
	return VAFString;
    }

    public double getPval() {
	if(pvalUncorr<-98) {
	    valueUnsetError("pval");
	}
	return pvalUncorr;
    }

    public double getMaxDiffVAF() {
	if(diffVAF<-98) {
	    valueUnsetError("diff");
	}
	return diffVAF;
    }

    public double getMaxVAF() {
	if(maxVAF<-98) {
	    valueUnsetError("max");
	}
	return maxVAF;
    }

    public int getNAmplicon() {
	if(namplicons<-98) {
	    valueUnsetError("n amplicons");
	}
	return namplicons;
    }

    public int getTotalBarcodes() {
	return bcmap.keySet().size();
    }

    public BaseFlagMap  aggregateOverBarcodes(PrintStream bcout) {

	boolean comma=false;

	BaseFlagMap overallMap=new BaseFlagMap();
	for(String rg : bcmap.keySet()) {
	    BaseFlagMap bfm=new BaseFlagMap();
	    for(String amp : bcmap.get(rg).keySet()) {
		bfm.add(bcmap.get(rg).get(amp));
	    }
            if(bfm.sum()>=minCountPerBC) {
	      char mc=bfm.maxBase();
	      int flag=bfm.maxFlag();
	      if(bfm.sum(mc)*1.0/bfm.sum()< minPercentRG ) {
		mc='N';
	      }
	      overallMap.add(new BaseFlag(mc, flag));
	      if(bcout!=null) {
		bfm.print();
		if(comma) {
		    bcout.print(",");
		}
		bcout.print(rg+"_"+mc);
		comma=true;
	      }
	  }
	}
	if(bcout!=null) {
	    bcout.println();
	}
	return overallMap;
    }
}




