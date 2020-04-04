import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class AminoAcidLLTester {


    @Test
    // Tests aminoAcidList
    public void aminoAcidList1(){
        char[] expected = {'R','K'};
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("AGGAAGCTGAGGGAGAGG");
        assertArrayEquals(expected, test.aminoAcidList());
    }
    @Test
    public void aminoAcidList2(){
        char[] expected = {'S'};
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("AGU");
        assertArrayEquals(expected, test.aminoAcidList());
    }
    @Test
    //test for amino counts
    public void aminoAcidCounts1() {
        int[] expected = {3};
        String testSequence = "GCGGCGGCGUGAAAGGU";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence(testSequence);
        assertArrayEquals(expected, test.aminoAcidCounts());

    }
    @Test
   // test for amino counts
    public void aminoAcidCounts2() {
        int[] expected = {1,2};
        String testSequence = "ACGGCGGCGUGAAAGGUAAAUUUUAGGGCC";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence(testSequence);
        assertArrayEquals(expected, test.aminoAcidCounts());

    }
    @Test
    // Tests createFromRNASequence to see if it works
    public void createFromRNASequence1(){
        String expected = "KGL";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("AAAGGCCUUUGAGAUAGAUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }

    @Test
    // Tests createFromRNASequence with stop codon in the middle
    public void createFromRNASequence2(){
        String expected = "";
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("CCC");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), test.aminoAcid);
            test = test.next;
        }
    }

    @Test
    //TEST FOR isSorted
    public void isSorted1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("AUUUUUGGGUAGUA");
        assertEquals(false, test.isSorted());
    }
    @Test
    public void isSorted2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GAAAAAUUUUUUGGGGAGGGGAUUGGGGUAAGGAUUUAG");
        assertEquals(false, test.isSorted());
    }
    @Test
    // Tests sort
    public void sort1(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GGG");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }
    @Test
    public void sort2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUUAAAGGGGUUUAAAGUA");
        test = AminoAcidLL.sort(test);
        assertEquals(true, test.isSorted());

    }

}