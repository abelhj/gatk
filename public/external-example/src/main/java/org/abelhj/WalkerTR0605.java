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

@Reference(window=@Window(start=-1, stop=1))
public class WalkerTR0605 extends LocusWalker<Integer,Integer>  {

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
    @Argument(fullName = "debug", shortName = "debug", doc= "1 for verbose output", required=false)
    int debug=0;
    ArrayList<Character> nucs=null;
    ArrayList<Integer> flags=null;
    
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
        //out.println("chr\tpos\tref\tdepth\trefct\tref83\tref163\tref99\tref147\tnumbc");
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
            String rg=p.getRead().getStringAttribute("X0");
            int flag=p.getRead().getFlags();
            char bb=(char)p.getBase();
            if(!bcreads.containsKey(rg)) {
                bcreads.put(rg, new ArrayList<BaseFlag>());
            }
            bcreads.get(rg).add(new BaseFlag(bb, flag));
            bfmap.add(bb, flag);           
        }
        
        char alt=bfmap.maxBase(refbase);
        if(alt=='N')
            alt='.';
        //out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+bfmap.sum(refbase)+"\t"+bfmap.printSums()+"\t");
        out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+bfmap.printSums(refbase)+"\t"+bfmap.printSums(alt)+"\t"+bfmap.printSums()+"\t");
        out.print(bcreads.keySet().size()+"\t");
	
        BaseFlagMap rgbfmap=new BaseFlagMap();

        for(String rg : bcreads.keySet()) {
            ArrayList<BaseFlag> reads=bcreads.get(rg);
            BaseFlagMap bfm=new BaseFlagMap();
            bfm.fill(reads);
            char mc=bfm.maxBase();
            int flag=reads.get(0).flag;
            if(!nucs.contains(mc) || bfm.sum(mc)<1 || bfm.sum(mc)/bfm.sum()<minPercentRG)
                mc='N';
            rgbfmap.add(mc, flag);
            if(debug==1)
                bfm.print();

        }
        char bcalt=rgbfmap.maxBase(refbase);
        if(bcalt=='N')
            bcalt='.';
        //out.print(alt+"\t"+rgbfmap.sum(refbase)+"\t"+rgbfmap.printSums()+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
        out.print(bcalt+"\t"+rgbfmap.printSums(refbase)+"\t"+rgbfmap.printSums(bcalt)+"\t"+rgbfmap.printSums()+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
       
        return 1;
    }
   


    /**
     * Provide an initial value for reduce computations. In this case we simply return an empty list
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;
    }

    /**
     * Outputs the number of genotypes called.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }

 
    /**
     * when we finish traversing, close the result list
     * @param result the final reduce result
     */
    public void onTraversalDone(Integer result) {
        
    }
}
