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

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;


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
    double minPercentRG=0.9;
    @Argument(fullName = "minCountPerBC", shortName = "minCtBC", doc= "Minimum number reads to accept bar code family", required=false)
    int minCountPerBC=1;
    @Argument(fullName = "debug", shortName = "debug", doc= "1 for verbose output", required=false)
    int debug=0;
    @Argument(fullName = "discardN", shortName = "discardN", doc= "alt is most common non-ref base excluding Ns", required=false)
    int discardN=0;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "writebc", shortName = "writebc", doc="0/1 to print list of barcodes", required=false)
    int writebc=0;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    boolean skipN=false;
    ArrayList<Character> nucs=null;
    ArrayList<Integer> flags=null;
    PrintStream bcout=null;
    
    public void initialize() {
	nucs=new ArrayList<Character>();
        nucs.add('A');
        nucs.add('C');
        nucs.add('G');
        nucs.add('T');
        nucs.add('N');
	flags=new ArrayList<Integer>();
        flags.add(83);
        flags.add(163);
        flags.add(99);
        flags.add(147);
	if(discardN==1) {
	    skipN=true;
	}
	if(writebc==1) {
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
	String loc=pileup.getLocation().toString();
        String [] chrpos=loc.split(":");
        LinkedHashMap<String, ArrayList<BaseFlag> > bcreads=new LinkedHashMap<String, ArrayList <BaseFlag> >();
        BaseFlagMap bfmap=new BaseFlagMap();
        
        
        Character refbase=(char)ref.getBase();
        
        for(PileupElement p : pileup) {
	    if(p.getRead().getIntegerAttribute("NM")<maxNM && p.getOffset()>=minOffset && p.getOffset()<=p.getRead().getReadLength()-minOffset) {
		String rg=p.getRead().getStringAttribute("X0");
		int flag=p.getRead().getFlags();
		char bb=(char)p.getBase();
		if(!bcreads.containsKey(rg)) {
		    bcreads.put(rg, new ArrayList<BaseFlag>());
		}
		bcreads.get(rg).add(new BaseFlag(bb, flag));
		bfmap.add(bb, flag);           
	    }
        }
        
        char alt=bfmap.maxBase(refbase, false);
        if(alt=='N') {
            alt='.';
	}
		
        out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+bfmap.printSums(refbase)+"\t"+bfmap.printSums(alt)+"\t"+bfmap.printSums()+"\t");
        out.print(bcreads.keySet().size()+"\t");
	if(writebc==1) {
	    bcout.print(chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t");
	}
	boolean comma=false;
	
        BaseFlagMap rgbfmap=new BaseFlagMap();

        for(String rg : bcreads.keySet()) {
            ArrayList<BaseFlag> reads=bcreads.get(rg);
            BaseFlagMap bfm=new BaseFlagMap();
            bfm.fill(reads);
            char mc=bfm.maxBase();
            int flag=reads.get(0).flag;
            if(!nucs.contains(mc) || bfm.sum(mc)<1 || bfm.sum(mc)/bfm.sum()<minPercentRG || bfm.sum()<minCountPerBC) {
                mc='N';
	    }
            rgbfmap.add(mc, flag);
            if(debug==1) {
                bfm.print();
	    }
	    if(writebc==1) {
		if(comma) {
		    bcout.print(",");
		}
		bcout.print(rg+"_"+mc);
		comma=true;
	    }
        }
	if(writebc==1) {
	    bcout.println();
	}
        char bcalt=rgbfmap.maxBase(refbase, skipN);
        if(bcalt=='N') {
            bcalt='.';
	}
	out.print(bcalt+"\t"+rgbfmap.printSums(refbase)+"\t"+rgbfmap.printSums(bcalt)+"\t"+rgbfmap.printSums()+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
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
