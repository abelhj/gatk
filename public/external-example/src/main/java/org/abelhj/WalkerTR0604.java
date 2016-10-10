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
public class WalkerTR0604 extends LocusWalker<Integer,Integer>  {

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
    
    public void initialize() {
        out.println("chr\tpos\tref\tdepth\trefct\tref83\tref163\tref99\tref147\tnumbc");
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
		String loc=pileup.getLocation().toString();
        String [] chrpos=loc.split(":");
	
		LinkedHashMap<Character, ArrayList<Character>> basect=new LinkedHashMap<Character, ArrayList<Character> >();
        LinkedHashMap<Character, LinkedHashMap<Integer,ArrayList<Integer> > > dirct=new LinkedHashMap<Character, LinkedHashMap<Integer, ArrayList<Integer> > >();
        LinkedHashMap<String, ArrayList<GATKSAMRecord> > bcreads=new LinkedHashMap<String, ArrayList <GATKSAMRecord> >();
        
        
        LinkedHashSet<Character> nucs=new LinkedHashSet<Character>();
        nucs.add('A');
        nucs.add('C');
        nucs.add('G');
        nucs.add('T');
        nucs.add('N');
        ArrayList<Integer> flags=new ArrayList<Integer>();
        flags.add(83);
        flags.add(163);
        flags.add(99);
        flags.add(147);
	
        for(Character c: nucs) {
            basect.put(c, new ArrayList<Character>() );
            LinkedHashMap<Integer, ArrayList<Integer> > temp=new LinkedHashMap<Integer, ArrayList<Integer> >();
            for(Integer fl: flags) {
                temp.put(fl, new ArrayList<Integer>());
            }
            dirct.put(c, temp);
        }
        
   
        
        for(PileupElement p : pileup) {
            String rg=p.getRead().getStringAttribute("X0");
            int flag=p.getRead().getFlags();
            if(!bcreads.containsKey(rg)) {
                bcreads.put(rg, new ArrayList<GATKSAMRecord>());
            }
            bcreads.get(rg).add(p.getRead());
            char bb=(char)p.getBase();
            if(nucs.contains(bb)) {
                basect.get(bb).add(bb);
                if(flags.contains(flag))
                    dirct.get(bb).get(flag).add(flag);
            }
        }
        Character refbase=(char)ref.getBase();
        out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+(char)ref.getBase()+"\t.\t"+pileup.getBases().length+"\t"+basect.get(refbase).size()+"\t");
        /*for(Integer fl:flags) {
            out.print(dirct.get(refbase).get(fl).size()+"\t");
        }*/
        for(Character c: nucs) {
            out.print(basect.get(c).size()+":");
            for(Integer fl:flags) 
                out.print(dirct.get(c).get(fl).size()+";");
            out.print("\t");
        }
        out.println(bcreads.keySet().size());
	return 1;
    }
        
        
      	
		/*if(totalnonref==0) {
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
       		
	    out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+(char)ref.getBase()+"\t"+bestalt+"\t"+bases.length+"\t"+totalnonref+"\t"+rgmap.keySet().size()+"\t"+rgct+"\t0\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");*/
   


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
