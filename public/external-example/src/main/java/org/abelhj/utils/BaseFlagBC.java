package org.abelhj.utils;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class BaseFlagBC extends BaseFlag{

    protected String barcode;
    protected int start;
    protected int end;

    public BaseFlagBC(char bb, int fl, String bc, int clipstart, int clipend) {
	super(bb, fl);
	barcode=bc;
	start=clipstart;
	end=clipend;
    }

    public BaseFlagBC(char bb, GATKSAMRecord rec) {
	super(bb, rec);
	barcode=rec.getStringAttribute("X0");
	start=rec.getSoftStart();
	end=rec.getSoftEnd();
    }

    public String toString() {
	String str="("+base+", "+flag+", "+barcode+", "+start+", "+end+")";
	return str;
    }

    public String getBarcode() {
	return barcode;
    }

    public int getStart() {
	return start;
    }

    public int getEnd() {
	return end;
    }
}