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
import org.abelhj.utils.ReadFamily;
import org.abelhj.utils.ReadFamilyLead;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;


import java.io.PrintStream;
import java.util.Arrays;
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
public class WalkerTRConsensus extends ReadWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    @Argument(fullName = "minPercentRG", shortName = "mprg", doc = "Minimum percent mutant bases per read group", required = false, minValue = 0, maxValue=1)
    double minPercentRG=0.9;
    @Argument(fullName = "debug", shortName = "debug", doc= "1 to print read groups and consensus to bcfile", required=false)
    int debug=0;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;
    @Argument(fullName = "consensusBam", shortName = "cbam", doc="consensus BAM ouput", required=true)
    String consensusBam=null;

    boolean debg=false;
    PrintStream bcout=null;
    SAMFileWriter samwriter=null;
    SAMFileWriterFactory sf=null;

    LinkedHashMap<String, LinkedHashMap<String, ReadFamilyLead> > bcmaster=null;
    LinkedHashMap<String, LinkedList<GATKSAMRecord> > bcreads;
    GenomeLoc oldpos=null;
    GenomeLoc curpos=null;
    
    public void initialize() {
  
	if(debug==1) {
	    debg=true;
	    if(bcfile==null) {
		bcfile=consensusBam+".debug.txt";
	    }
	    try {
		bcout=new PrintStream(new File(bcfile));
	    } catch(Exception e) {
		System.err.println("Could be creat debug output file.\n");
	    }
	}
	bcmaster=new LinkedHashMap<String, LinkedHashMap<String, ReadFamilyLead> >();
	bcreads=new LinkedHashMap<String, LinkedList<GATKSAMRecord> >();
	sf=new SAMFileWriterFactory();
	//samwriter=ReadUtils.createSAMFileWriterWithCompression(this.getToolkit().getSAMFileHeader(), true, consensusBam, 5);
	samwriter=sf.makeBAMWriter(this.getToolkit().getSAMFileHeader(), true, new File(consensusBam), 5);
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
	
	if(ref!=null) {
	    if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case
	    curpos=ref.getLocus();
	    if(!curpos.equals(oldpos) ) {
		for(String bc : bcreads.keySet()) {
		    ReadFamily rf=new ReadFamily(bcreads.get(bc));
		    rf.getConsensus(bcout, bcmaster, samwriter);	
		}
		bcreads=new LinkedHashMap<String, LinkedList<GATKSAMRecord> >();
		oldpos=curpos;
	    }
	    String bc=read.getStringAttribute("X0");
	    if(!bcreads.containsKey(bc)) {
		bcreads.put(bc, new LinkedList<GATKSAMRecord>());
	    }
	    bcreads.get(bc).add(read);
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
	for(String bc : bcreads.keySet()) {
	    ReadFamily rf=new ReadFamily(bcreads.get(bc));
	    rf.getConsensus(bcout, bcmaster, samwriter);
	}
	samwriter.close();
    }
}
