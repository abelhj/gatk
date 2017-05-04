package org.abelhj;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;


import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.PrintStream;

import org.abelhj.utils.BaseFlag;
import org.abelhj.utils.BaseFlagMap;
import org.abelhj.utils.BaseFlagBC;
import org.abelhj.utils.BaseFlagMapBCNoAmp;


@Reference(window=@Window(start=-1, stop=1))
public class WalkerTR1203 extends LocusWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    @Argument(fullName = "minPercentRG", shortName = "mprg", doc = "Minimum percent mutant bases per read group", required = false, minValue = 0, maxValue=1)
    double minPercentRG = 0.9;
    @Argument(fullName = "minCountPerBC", shortName = "minCtBC", doc= "Minimum number reads to accept bar code family", required=false)
    int minCountPerBC = 1;
    @Argument(fullName = "debug", shortName = "debug", doc= "true for verbose output", required=false)
    boolean debug = false;
    @Argument(fullName = "allowN", shortName = "allowN", doc= "allow N as alt", required=false)
    boolean allowN = false;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile = null;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset = 0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    PrintStream bcout=null;
    
    public void initialize() {

	if(debug) {
	    if(bcfile!=null) {
		try {
		    bcout=new PrintStream(new File(bcfile));
		} catch(Exception e) {
		    System.err.println("barcode file not found\n");
		}
	    } else {
		System.err.println("Error: Must provide name for barcode list file.\n");
		System.exit(1);
	    }
	}
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	GenomeLoc loc=pileup.getLocation();
        char refbase=(char)ref.getBase();
	BaseFlagMap bfmap=new BaseFlagMap();
	BaseFlagMapBCNoAmp ecmap=new BaseFlagMapBCNoAmp(loc, refbase, minPercentRG, minCountPerBC);
        
        for(PileupElement p : pileup) {

	    GATKSAMRecord pread=p.getRead();
	    if(p.getRead().getReadPairedFlag() && p.getRead().getProperPairFlag()) {
		if(pread.getIntegerAttribute("NM")<maxNM && p.getOffset()>=minOffset && p.getOffset()<=pread.getReadLength()-minOffset) {
		    BaseFlagBC bfl=new BaseFlagBC((char)p.getBase(), pread);
		    ecmap.add(bfl);
		    bfmap.add(bfl);           
		}
	    }
	}

	char alt=bfmap.maxBase(refbase, false);
	char altEC='.';
	double vaf=0;
	double vafEC=0;
	BaseFlagMap bfmapEC = new BaseFlagMap();

	if(alt=='N') {
	    alt='.';
	    bfmapEC=new BaseFlagMap();
	} else {
	    vaf=bfmap.calcVAF(refbase, alt);
	    bfmapEC=ecmap.aggregateOverBarcodes(bcout);
	    altEC=bfmapEC.maxBase(refbase, !allowN);
	    if(altEC=='N') {
		altEC='.';
	    } else {
		vafEC=bfmapEC.calcVAF(refbase, altEC);
	    }
	}
	String [] chrpos=loc.toString().split(":");
	String nonbcstr=chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+String.format("%.4e", vaf)+"\t"+bfmap.printSums(refbase)+"\t"+bfmap.printSums(alt)+"\t"+bfmap.printSums();
	String bcstr=altEC+"\t"+bfmapEC.sum()+"\t"+String.format("%.4e", vafEC)+"\t"+bfmapEC.printSums(refbase)+"\t"+bfmapEC.printSums(altEC)+"\t"+bfmapEC.printSums();
	System.out.print(nonbcstr+"\t"+bcstr+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
        return 1;
    }
   
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }

    public void onTraversalDone(Integer result) {
        
    }
}
