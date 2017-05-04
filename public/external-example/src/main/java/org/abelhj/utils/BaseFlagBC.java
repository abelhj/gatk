package org.abelhj.utils;

public class BaseFlagBC extends BaseFlag{

    protected String barcode;
    protected int start;
    protected int end;

    public BaseFlagBC(char bb, int fl, String bc, int clipstart , int clipend) {
	super(bb, fl);
	barcode=bc;
	start=clipstart;
	end=clipend;
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