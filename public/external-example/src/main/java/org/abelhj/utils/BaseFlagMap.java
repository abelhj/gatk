package org.abelhj.utils;

import java.util.Map;
import java.util.Array;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import htsjdk.samtools.SAMFlag;

public class BaseFlagMap {


    public List<Character> nucs=null;
    public List<Integer> flags=null;

    private Map<Character, Map<Integer, Integer > > map=null;


    public BaseFlagMap(){
        nucs=new ArrayList<Character>(Arrays.asList('A', 'C', 'G', 'T', 'N'));
        flags=new ArrayList<Integer>(Arrays.asList(
						   SAMFlag.FIRST_OF_PAIR.intValue(),
						   SAMFlag.SECOND_OF_PAIR.intValue(),
						   SAMFlag.FIRST_OF_PAIR.intValue()+SAMFlag.READ_REVERSE_STRAND.intValue(),
						   SAMFlag.SECOND_OF_PAIR.intValue()+SAMFlag.READ_REVERSE_STRAND.intValue()));
        map=new LinkedHashMap<Character, Map<Integer,Integer > >();
        for(Character c : nucs) {
            Map<Integer, Integer > temp=new LinkedHashMap<Integer, Integer >();
            for(Integer fl : flags) {
                temp.put(fl, 0);
            }
            map.put(c, temp);
        }
    }

    public void fill(List<BaseFlag> list) {
	for(BaseFlag val : list) 
	    add(val);
    }

    public boolean add(BaseFlag val) {
        if(nucs.contains(val.getBase()) && flags.contains(val.getFlag())) {
            map.get(val.getBase()).put(val.getFlag(), map.get(val.getBase()).get(val.getFlag())+1);
            return true;
        } else return false;
    }


    public void add(BaseFlagMap bfm) {
	for (Character c : nucs) {
	    for (Integer f : flags ) {
		int oldval=map.get(c).get(f);
		map.get(c).put(f, oldval+bfm.sum(c, f));
	    }
	}
    }


    public int sum() {
	int ss=0;
	for(Character c : nucs) {
	    for(Integer f : flags) 
		ss+=map.get(c).get(f);
	}
	return ss;
    }
	
    public int sum(char base) {
        if(nucs.contains(base)) {
            int ss=0;
            for(Integer fl:flags) 
                ss+=map.get(base).get(fl);
            return ss;
        }
        else {
            System.err.println("warning: unexpected base "+base);
            return 0;
        }
    }

    public int sum(char base, int flag) {
        if(nucs.contains(base) && flags.contains(flag)) {
            return map.get(base).get(flag);
        } else {
            System.err.println("Warning: unexpected base/flag combination "+base+" "+flag);
            return 0;
        }
    }

    public double calcVAF(char ref, char alt) {
	int refct=sum(ref);
	int altct=sum(alt);
	double vaf=1.0*altct/(refct+altct);
	return vaf;
    }

    public void print() {
	for(Character c: nucs) {
	    for(Integer fl : flags) {
		System.err.println(c+"\t"+fl+"\t"+map.get(c).get(fl));
	    }
	}
    }
    
    public String printSums() {
        String ret="";
        for(Character c : nucs) {
            ret+=printSums(c);
        }
        return ret;
    }

    public String printSums(char c) {
        String ret="";
        if(nucs.contains(c)) {
            ret+=sum(c)+":";
	    int or1=sum(c, SAMFlag.FIRST_OF_PAIR.intValue())+sum(c, SAMFlag.SECOND_OF_PAIR.intValue()+SAMFlag.READ_REVERSE_STRAND.intValue());
            int or2=sum(c, SAMFlag.SECOND_OF_PAIR.intValue())+sum(c,  SAMFlag.FIRST_OF_PAIR.intValue()+SAMFlag.READ_REVERSE_STRAND.intValue());
            ret+=or1+","+or2+";";
        } else {
            ret+=".:.,.;";
        }
        return ret;
    }

    public int maxFlag() {
	int maxct=0;
	int maxflag=0;
	for(Integer fl : flags) {
	    int ct=0;
	    for(Character c : nucs ) {
		ct+=map.get(c).get(fl);
	    }
	    if(ct>maxct) {
		maxct=ct;
		maxflag=fl;
	    }
	}
	return maxflag;
    }

    public char maxBase() {
	int maxct=0;
	char maxbase='N';
	for(Character c : nucs) {
	    int ss=sum(c);
	    if(ss>maxct) {
		maxct=ss;
		maxbase=c;
	    }
	}
	return maxbase;
    }

    public char maxBase(boolean ignoreN) {
	if(!ignoreN) {
	    return maxBase();
	} else {
	  int maxct=0;
	  char maxbase='N';
	  for(Character c : nucs) {
	      if(c!='N') {
		  int ss=sum(c);
		  if(ss>maxct) {
		      maxct=ss;
		      maxbase=c;
		  }
	      }
	  }
	  return maxbase;
	}
    }

    public char maxBase(char ref) {
	int maxct=0;
	char maxbase='N';
	for(Character c : nucs) {
	    if(c!=ref) {
		int ss=sum(c);
		if(ss>maxct) {
		    maxct=ss;
		    maxbase=c;
		}
	    }
	}
        return maxbase;
    }
    
    public char maxBase(char ref, boolean ignoreN) {
	if(!ignoreN) {
	    return maxBase(ref);
	} else {
	    int maxct=0;
	    char maxbase='N';
	    for(Character c : nucs) {
		if(c!=ref && c!='N') {
		    int ss=sum(c);
		    if(ss>maxct) {
			maxct=ss;
			maxbase=c;
		    }
		}
	    }
	    return maxbase;
	}
    }
}




