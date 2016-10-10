package org.abelhj.utils;

public class BaseFlag1{

    public char base;
    public int flag;

    public BaseFlag1(char bb, int fl) {
	base=bb;
	flag=fl;
    }

    public String toString() {
	String str="("+base+", "+flag+")";
	return str;
    }
}