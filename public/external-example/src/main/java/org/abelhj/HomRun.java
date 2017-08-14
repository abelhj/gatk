package org.abelhj;
//package org.broadinstitute.gatk.utils.examples;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.abelhj.utils.TypedTuple;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileWriterFactory;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;


@Reference(window=@Window(start=-1, stop=1))
public class HomRun extends ReadWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */

    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    //@Argument(fullName = "debug", shortName = "debug", doc= "1 to print read groups and consensus to bcfile", required=false)
    //boolean debug=false;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=12;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    int pos=31022442;

    int[] counts=new int[20];
    
    public void initialize() {
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
	
	if(ref!=null) {
            if(read.getFirstOfPairFlag() && read.getIntegerAttribute("NM")<maxNM && read.getAlignmentStart()<pos && read.getAlignmentEnd()>pos ) {
		int readpos=ReadUtils.getReadCoordinateForReferenceCoordinate(read, pos, ReadUtils.ClippingTail.LEFT_TAIL);
		if(readpos>minOffset && readpos < read.getReadLength()-minOffset) {
		    int len=getRunLength(readpos, read);
		    if(len<20) {
			counts[len]++;
		    }
		}
	    }
	}
	return 1;
    }
   
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }

    public void onTraversalDone(Integer result) {
	for(int ii=0; ii<counts.length; ii++) {
	    System.out.println(ii+"\t"+counts[ii]);
	}
    }


    public static int getRunLength(int readpos, GATKSAMRecord read) {

	int maxlen=10;
	char targetbase='G';
	byte[] bases=read.getReadBases();
	int offset=readpos;

	int len=0;
	int left=offset;
	int right=offset;
	if((char) bases[offset]==targetbase) {
	    while((char)bases[left-1]==targetbase && left>offset-maxlen) {
		left--;
	    }
	    while((char)bases[right+1]==targetbase && right<offset+maxlen) {
		right++;
	    }
	    len=right-left+1;
	}
	System.err.println(read.getReadString()+"\t"+offset+"\t"+left+"\t"+right+"\t"+len);
	return len;
    }
}
