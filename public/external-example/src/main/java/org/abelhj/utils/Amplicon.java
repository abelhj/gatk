package org.abelhj.utils;

import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class Amplicon {
   
    protected int start;
    protected int end;
    protected String name;
    protected String chr;
    protected char strand;

    public Amplicon(String chr, int start, int end, String nm, char strand) {
        this.chr=chr;
	this.start=start;
	this.end=end;
        this.name=nm;
	this.strand=strand;
    }

    public Amplicon(String chr, int start, int end, String nm) {
	this(chr, start, end, nm, '+');
    }

    public String toString() {
	String str=chr+":"+start+"-"+end+":"+name+":"+strand;
        return str;
    }

    public int getStart() {
	return start;
    }

    public int getEnd() {
	return end;
    }

    public String getChrom() {
	return chr;
    }

    public String getName() {
	return name;
    }

    public char getStrand() {
	return strand;
    }
	

    public int getDist(GATKSAMRecord rec) {
	int dist=0;
	if(rec.getReferenceName().equals(chr)) {
	    dist=Math.abs(rec.getAlignmentStart()-start);
            int dist2=Math.abs(rec.getMateAlignmentStart()-start);
	    if(dist2<dist)
		dist=dist2;
	    return dist;
	} else {
	    return 9999;
	}
    }


}