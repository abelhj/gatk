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
import java.util.List;

import org.abelhj.utils.BaseFlag;
import org.abelhj.utils.BaseFlagMap;
import org.abelhj.utils.BaseFlagBC;
import org.abelhj.utils.BaseFlagMapBC;


@Reference(window=@Window(start=-1, stop=1))
public class HomoLen extends LocusWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=12;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    PrintStream bcout=null;
    
    public void initialize() {

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

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
	
	int maxlen=12;
        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	GenomeLoc loc=pileup.getLocation();
        char refbase=(char)ref.getBase();
        char targetbase='G';
	List<Integer> counts=new ArrayList<Integer>();
        for(PileupElement p : pileup) {

	    GATKSAMRecord pread=p.getRead();
            int offset=p.getOffset();  //6 to 11
	    if(pread.getIntegerAttribute("NM")<maxNM && offset>=minOffset && offset<=pread.getReadLength()-minOffset) {
		int len=0;
		int left=offset;
		int right=offset;
                if((char)p.getBase()==targetbase) {
                    byte[] bases=pread.getReadBases();
		    while((char)bases[left]==targetbase && left>offset-maxlen) {
			left--;
		    }
		    while((char)bases[right]==targetbase && right<offset+maxlen) {
			right++;
		    }
		    len=right-left+1;
		}
		System.err.println(pread.getReadString()+"\t"+offset+"\t"+left+"\t"+right+"\t"+len);
		counts.set(len, counts.get(len)+1);
	    }
        }
	for(int i=0; i<counts.size(); i++) {
	    System.out.println(i+"\t"+counts.get(i));
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
    }
}
