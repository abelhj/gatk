package org.abelhj;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;

import htsjdk.samtools.SAMFlag;

import java.io.PrintStream;
import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;

import org.abelhj.utils.BaseFlagLocMap;
import org.abelhj.utils.BaseFlagMap;
import org.abelhj.utils.BaseFlagLoc;
import org.abelhj.utils.BaseFlag;

@Reference(window=@Window(start=-1, stop=1))
public class WalkerTR2017_1 extends LocusWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    @Argument(fullName = "minPercentRG", shortName = "mprg", doc = "Minimum percent mutant bases per read group", required = false, minValue = 0, maxValue=1)
    double minPercentRG=0.9;
    @Argument(fullName = "minCountPerBC", shortName = "minCtBC", doc= "Minimum number reads to accept bar code family", required=false)
    int minCountPerBC=1;
    @Argument(fullName = "debug", shortName = "debug", doc= "1 for verbose output", required=false)
    int debug=0;
    @Argument(fullName = "discardN", shortName = "discardN", doc= "alt is most common non-ref base excluding Ns", required=false)
    int discardN=0;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    boolean skipN=false;
    PrintStream bcout=null;
    boolean debugBool=false;
    int minPerAmp=50;
    
    public void initialize() {
	if(discardN==1) {
	    skipN=true;
	}
	if(debug==1) {
	    debugBool=true;
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
	BaseFlagLocMap ecmap=new BaseFlagLocMap( loc, refbase, minPerAmp, minPercentRG, minCountPerBC);
        
        for(PileupElement p : pileup) {

	    GATKSAMRecord pread=p.getRead();
	    if(pread.getIntegerAttribute("NM")<maxNM && p.getOffset()>=minOffset && p.getOffset()<=pread.getReadLength()-minOffset) {

		String rg=pread.getStringAttribute("X0");
		int flag=0;
		if(pread.getFirstOfPairFlag()) {
		    flag+=SAMFlag.FIRST_OF_PAIR.intValue();
		}
		if(pread.getSecondOfPairFlag()) {
		    flag+=SAMFlag.SECOND_OF_PAIR.intValue();
		}
		if(pread.getReadNegativeStrandFlag()) {
		    flag+=SAMFlag.READ_REVERSE_STRAND.intValue();
		} 
		BaseFlagLoc bfl=new BaseFlagLoc((char)p.getBase(), flag, pread.getSoftStart(), pread.getSoftEnd());
		ecmap.add(bfl, rg);
		bfmap.add(bfl);           
	    }
        }

	char alt=bfmap.maxBase(refbase, false);
	char altEC='.';
	double vaf=0;
	double vafEC=0;
	BaseFlagMap bfmapEC=null;
	boolean ampBias=false;

	double pval=1;
	double maxDiff=0;
	int namplicons=0;

	if(alt=='N') {
	    alt='.';
	    bfmapEC=new BaseFlagMap();
	} else {
	    vaf=bfmap.calcVAF(refbase, alt);
	    ecmap.setAlt(alt);
	    ecmap.calcAmpBias(debugBool);
	    pval=ecmap.getPval();
	    namplicons=ecmap.getNAmplicon();
	    maxDiff=ecmap.getMaxDiffVAF();
	    if(pval<1e-5 && maxDiff>0.05) {
		ampBias=true;
	    }
	    bfmapEC=ecmap.calcOverallVAF(ampBias, bcout, debugBool);
	    altEC=bfmapEC.maxBase(refbase, false);
	    if(altEC=='N') {
		altEC='.';
	    } else {
		vafEC=bfmapEC.calcVAF(refbase, altEC);
	    }
	}
	String [] chrpos=loc.toString().split(":");
	String nonbcstr=chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+String.format("%.4e", vaf)+"\t"+bfmap.printSums(refbase)+"\t"+bfmap.printSums(alt)+"\t"+bfmap.printSums();
	String bcstr=altEC+"\t"+bfmapEC.sum()+"\t"+String.format("%.4e", vafEC)+"\t"+bfmapEC.printSums(refbase)+"\t"+bfmapEC.printSums(altEC)+"\t"+bfmapEC.printSums();
	System.out.print(nonbcstr+"\t"+bcstr+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\t"+ampBias+"\t"+maxDiff+"\t"+pval+"\t"+namplicons+"\n");
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
