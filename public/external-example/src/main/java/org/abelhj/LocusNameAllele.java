package org.abelhj;

import java.util.LinkedHashMap;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

public class LocusNameAllele {
    private GenomeLoc locus=null;
    private LinkedHashMap<String, PileupElement> reads=null;

    public LocusNameAllele(GenomeLoc loc, LinkedHashMap<String, PileupElement> map) {
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
    
    public LinkedHashMap<String, PileupElement> getAlleles() {
        return reads;
    }
    
    
}