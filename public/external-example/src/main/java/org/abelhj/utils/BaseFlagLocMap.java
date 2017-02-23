package org.abelhj.utils;

import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.PrintStream;

import org.broadinstitute.gatk.utils.GenomeLoc;

import htsjdk.samtools.SAMFlag;

public class BaseFlagLocMap {

    private LinkedHashMap<Integer, Integer > endPos=null;
    private LinkedHashMap<Integer, LinkedHashMap<String, BaseFlagMap> > map=null;
    private LinkedHashMap<String, LinkedHashMap<Integer, BaseFlagMap> > bcmap=null;
    private int minBCPerAmp=50;
    private char refBase='\0';
    private char altBase='\0';
    private double minPercentRG=0.8;
    private int minCountPerBC=1;
    private GenomeLoc loc=null;
    private double diffVAF=0;
    private double pvalUncorr=1;
    private int refPos;
    private int namplicons=0;

    public BaseFlagLocMap(GenomeLoc pos, char ref, int minPerAmp, double minPct, int minCount){
	refBase=ref;
	map=new LinkedHashMap<Integer, LinkedHashMap<String, BaseFlagMap> >();
	bcmap=new LinkedHashMap<String, LinkedHashMap<Integer, BaseFlagMap> >();
	endPos=new LinkedHashMap<Integer, Integer>();
	loc=pos;
	refPos=loc.getStart();
	minBCPerAmp=minPerAmp;
	minPercentRG=minPct;
	minCountPerBC=minCount;
    }
    
    public void setAlt(char alt) {
	altBase=alt;
    }

    public void add(BaseFlagLoc bfl, String barcode) {

	BaseFlag bf=new BaseFlag(bfl.base, bfl.flag);
	if(map.containsKey(bfl.start) && map.get(bfl.start).containsKey(barcode)) {
	    map.get(bfl.start).get(barcode).add(bf);
	} else {
	    BaseFlagMap bfm=new BaseFlagMap();
	    bfm.add(bf);

	    if(!map.containsKey(bfl.start)) {
		map.put(bfl.start, new LinkedHashMap<String, BaseFlagMap>());
		endPos.put(bfl.start, bfl.end);
	    } 
	    if(!map.get(bfl.start).containsKey(barcode)) {
		map.get(bfl.start).put(barcode, bfm);
	    }

	    if(!bcmap.containsKey(barcode)) {
		bcmap.put(barcode, new LinkedHashMap<Integer, BaseFlagMap>());
	    }
	    if(!bcmap.get(barcode).containsKey(bfl.start)) {
		bcmap.get(barcode).put(bfl.start, bfm);
	    }
	}
    }


    public void calcAmpBias(boolean debug) {

	ArrayList<Integer> refct=new ArrayList<Integer>();
	ArrayList<Integer> altct=new ArrayList<Integer>();
	ArrayList<Integer> ends=new ArrayList<Integer>();
	
	for (int cstart : map.keySet()) {
	    if(map.get(cstart).keySet().size()>=minBCPerAmp) {
	       	//System.err.println(cstart);
		BaseFlagMap bf=new BaseFlagMap();
		for(String bc : map.get(cstart).keySet()) {
		    
		    BaseFlagMap bfm=map.get(cstart).get(bc);
		    char maxBase=map.get(cstart).get(bc).maxBase( false);
		    int maxCt=bfm.sum(maxBase);
		    int total=bfm.sum();

		    if( maxCt*1.0/total > minPercentRG && total > minCountPerBC) {
			bf.add(new BaseFlag(maxBase, SAMFlag.FIRST_OF_PAIR.intValue()));
		    }		
		}
		refct.add(bf.sum(refBase));
		altct.add(bf.sum(altBase));
		ends.add(endPos.get(cstart));
	    }
	}
	diffVAF=ChiSquareUtils.maxDiffVAF(refct, altct);
	pvalUncorr=ChiSquareUtils.chiSquare(refct, altct);
	namplicons=refct.size();
	if(debug) {
	  System.err.println(diffVAF+"\t"+pvalUncorr);
	  for(int i=0; i<refct.size(); i++) {
	      System.err.println(i+"\t"+refct.get(i)+"\t"+altct.get(i)+"\t"+ends.get(i));
	  }
	}
    }

    public double getPval() {
	return pvalUncorr;
    }

    public double getMaxDiffVAF() {
	return diffVAF;
    }

    public int getNAmplicon() {
	return namplicons;
    }

    public BaseFlagMap  calcOverallVAF(boolean biased, PrintStream bcout, boolean debug) {

	boolean comma=false;
	int minDist=20;

	BaseFlagMap overallMap=new BaseFlagMap();
	for(String rg : bcmap.keySet()) {
	    BaseFlagMap bfm=new BaseFlagMap();
	    // if there is amplicon bias, ignore readgroup where variant is too close to either soft-clipped end.
	    if(biased) {     
		int start=bcmap.get(rg).keySet().iterator().next();
		int end=endPos.get(start);
		if(Math.abs(start-refPos)<minDist || Math.abs(end-refPos)<minDist) {
		    continue;
		}
	    }
	    for(int cstart : bcmap.get(rg).keySet()) {
		bfm.add(bcmap.get(rg).get(cstart));
	    }
            if(bfm.sum()>=minCountPerBC) {
	      char mc=bfm.maxBase();
	      int flag=bfm.maxFlag();
	      if(bfm.sum(mc)*1.0/bfm.sum()< minPercentRG ) {
		mc='N';
	      }
	      overallMap.add(new BaseFlag(mc, flag));
	      if(debug) {
		bfm.print();
		if(comma) {
		    bcout.print(",");
		}
		bcout.print(rg+"_"+mc);
		comma=true;
	      }
	  }
	}
	if(debug) {
	    bcout.println();
	}
	return overallMap;
    }
}




