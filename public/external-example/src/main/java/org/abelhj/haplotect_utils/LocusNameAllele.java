package org.abelhj.haplotect_utils;

import java.util.LinkedHashMap;

import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.GenomeLoc;

public class LocusNameAllele {
    private GenomeLoc locus=null;
    private LinkedHashMap<String, BaseandQual> reads=null;

    public LocusNameAllele(GenomeLoc loc, LinkedHashMap<String, BaseandQual> map) {
        locus=loc;
        reads=map;
    }
    
    public String toString() {
        String ret=locus+"\t";
        for(String str:reads.keySet()) {
            ret+=str+":"+reads.get(str).getBase()+"\n\t\t";
        }
        return ret;
    }
    
    public GenomeLoc getLocus() {
        return locus;
    }
    
    public LinkedHashMap<String, BaseandQual> getAlleles() {
        return reads;
    }
    
    
}