package org.abelhj.utils;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class BaseFlag {

    protected char base;
    protected int flag;

    public BaseFlag(char bb, int fl) {
	base=bb;
	flag=fl;
    }

    public BaseFlag(char bb, GATKSAMRecord rec) {
	base=bb;
	flag=0;
	if(rec.getFirstOfPairFlag()) {
	    flag+=SAMFlag.FIRST_OF_PAIR.intValue();
	}
	if(rec.getSecondOfPairFlag()) {
	    flag+=SAMFlag.SECOND_OF_PAIR.intValue();
	}
	if(rec.getReadNegativeStrandFlag()) {
	    flag+=SAMFlag.READ_REVERSE_STRAND.intValue();
	}
    }	

    public String toString() {
	String str="("+base+", "+flag+")";
	return str;
    }

    public char getBase() {
	return base;
    }

    public int getFlag() {
	return flag;
    }
}