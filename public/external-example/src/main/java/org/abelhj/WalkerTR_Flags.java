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

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;

@Reference(window=@Window(start=-1, stop=1))
public class WalkerTR_Flags extends LocusWalker<Integer,Integer>  {

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

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

	ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	byte bases[] = pileup.getBases();
	int totalnonref=0;
	for(int index=0; index<bases.length; index++) {
	    if(bases[index]!=ref.getBase())
		totalnonref++;
	}
	String loc=pileup.getLocation().toString();
	String [] chrpos=loc.split(":");
	
	LinkedHashMap<String , int[]> rgmap=new LinkedHashMap<String, int[]>();
	LinkedHashMap<String, ArrayList<Integer> > rgflags= new LinkedHashMap<String, ArrayList<Integer> > ();
	LinkedHashMap<Character, Integer> altct=new LinkedHashMap<Character, Integer>();
	LinkedHashSet<Character> nucs=new LinkedHashSet<Character>();
	nucs.add('A');
	nucs.add('C');
	nucs.add('G');
	nucs.add('T');
	nucs.add('N');
	
	for(Character c: nucs) 
		altct.put(c, 0);
      	
	for(PileupElement p : pileup) {
		String rg=p.getRead().getStringAttribute("X0");
		//System.err.print(p.getRead().getFlags());
		if(!rgmap.containsKey(rg)) {
			rgmap.put(rg, new int[2]);
			rgflags.put(rg, new ArrayList<Integer>());
		}
		rgflags.get(rg).add(p.getRead().getFlags());
		if(p.getBase() == ref.getBase() )
			rgmap.get(rg)[0]++;
		else {
			rgmap.get(rg)[1]++;
			altct.put((char)p.getBase(), altct.get((char)p.getBase())+1);
		}
	}
	for(String s: rgflags.keySet()) {
	    System.err.print(s+"\t");
	    for(Integer k:rgflags.get(s)){
		System.err.print(k+" ");
	    }
	    System.err.println();
	}
	System.err.println();
	if(totalnonref==0) {
	    out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+(char)ref.getBase()+"\t.\t"+bases.length+"\t"+totalnonref+"\t"+rgmap.keySet().size()+"\t0\t0\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
      	} else {

	    int rgct=0;
	    for(String rg: rgmap.keySet()) {
		if((double)rgmap.get(rg)[1]/(rgmap.get(rg)[0]+rgmap.get(rg)[1])>minPercentRG) 
		    rgct++;
	    }
	    char bestalt='N';
	    int totalalt=0;
	    for(Character c: nucs) 
		totalalt+=altct.get(c);
	    if(totalalt>0) {
		for(Character c: nucs) {
		    if((double)altct.get(c)/totalalt>0.9)
				bestalt=c;
		}
	    }
       		
	    out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+(char)ref.getBase()+"\t"+bestalt+"\t"+bases.length+"\t"+totalnonref+"\t"+rgmap.keySet().size()+"\t"+rgct+"\t0\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
	}
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
