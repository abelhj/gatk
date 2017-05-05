package org.abelhj.utils;

import java.util.List;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.PrintStream;

import org.broadinstitute.gatk.utils.GenomeLoc;

public class BaseFlagMapBCNoAmp {

    private Map<String, BaseFlagMap> bcmap=null;
    private char refBase='\0';
    private double minPercentRG=0.8;
    private int minCountPerBC=1;
    private GenomeLoc loc=null;

    public BaseFlagMapBCNoAmp(GenomeLoc pos, char ref, double minPct, int minCount){
	refBase=ref;
	bcmap=new LinkedHashMap<String,  BaseFlagMap >();
	loc=pos;
	minPercentRG=minPct;
	minCountPerBC=minCount;
    }
    
    public void add(BaseFlagBC bfl) {

	if(bcmap.containsKey(bfl.getBarcode())) {
	    bcmap.get(bfl.getBarcode()).add(bfl);
	} else {
	    BaseFlagMap bfm=new BaseFlagMap();
	    bfm.add(bfl);
	    bcmap.put(bfl.getBarcode(), bfm);
	}
    }

    public BaseFlagMap aggregateOverBarcodes ( PrintStream bcout) {

	boolean comma=false;

	BaseFlagMap overallMap=new BaseFlagMap();
	for(String rg : bcmap.keySet()) {
	    BaseFlagMap bfm=bcmap.get(rg);
            if( bfm.sum()>=minCountPerBC ) {
	      char mc=bfm.maxBase();
	      int flag=bfm.maxFlag();
	      if(bfm.sum(mc)*1.0/bfm.sum()< minPercentRG || bfm.sum(mc)<1 ) {
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




