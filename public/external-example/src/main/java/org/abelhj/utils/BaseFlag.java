package org.abelhj.utils;

public class BaseFlag {

    protected char base;
    protected int flag;

    public BaseFlag(char bb, int fl) {
	base=bb;
	flag=fl;
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