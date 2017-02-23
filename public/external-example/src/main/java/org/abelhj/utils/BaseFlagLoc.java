package org.abelhj.utils;

public class BaseFlagLoc{

    public char base;
    public int flag;
    public int start;
    public int end;

    public BaseFlagLoc(char bb, int fl, int clipstart , int clipend) {
	base=bb;
	flag=fl;
	start=clipstart;
	end=clipend;
    }

    public String toString() {
	String str="("+base+", "+flag+", "+start+", "+end+")";
	return str;
    }
}