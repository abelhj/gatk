package org.abelhj.utils;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class BaseFlagAmpBC extends BaseFlag{

    protected String barcode;
    protected String ampName;

    public BaseFlagAmpBC(char bb, int fl, String bc, String name) {
	super(bb, fl);
	barcode=bc;
	ampName=name;
    }

    public BaseFlagAmpBC(char bb, GATKSAMRecord rec) {
	super(bb, rec);
	barcode=rec.getStringAttribute("X0");
	ampName=rec.getStringAttribute("X1");
    }

    public String toString() {
	String str="("+base+", "+flag+", "+barcode+", "+ampName+")";
	return str;
    }

    public String getBarcode() {
	return barcode;
    }

    public String getAmp() {
	return ampName;
    }
}