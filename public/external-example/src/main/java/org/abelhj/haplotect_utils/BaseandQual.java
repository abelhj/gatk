package org.abelhj.haplotect_utils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

public class BaseandQual{

    private byte base;
    private byte qual;
   
    public BaseandQual(PileupElement pel) {
	base=pel.getBase();
	qual=pel.getQual();
    }

    public byte getBase() {
	return base;
    }

    public byte getQual() {
	return qual;
    }
}