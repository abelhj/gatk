package org.abelhj.utils;

public class BaseFlag {

    public char base;
    public int flag;

    public BaseFlag(char bb, int fl) {
	base=bb;
	flag=fl;
    }

    public String toString() {
	String str="("+base+", "+flag+")";
	return str;
    }
}