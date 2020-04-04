class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;


  AminoAcidLL() {
  }

  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon) {

    this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    this.codons = AminoAcidResources.getCodonListForAminoAcid(this.aminoAcid);
    this.counts = new int[codons.length];
    incC(inCodon);
    this.next = null;
  }

  private void incC(String inCodon) {
    for (int i = 0; i < this.codons.length; i++) {
      if (this.codons[i].equals(inCodon)) {
        this.counts[i]++;
      }
    }
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops,
   * if not passes the task to the next node.
   * If there is no next node, add a new node to the list that would contain the codon.
   */
  private void addCodon(String inCodon) {
    //compares the current amino acid with results
    if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
      this.incC(inCodon);
    }
    //recursively goes to the next
    else if (next != null) {
      next.addCodon(inCodon);
    }
    //In the end  we make a new node
    else {
      this.next = new AminoAcidLL(inCodon);
    }
  }




  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for (int i = 0; i < this.counts.length; i++) {
      sum += counts[i];
    }
    return sum;
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts.
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList) {
   //base case
    if (this.next == null && inList.next == null){
      return (this.aminoAcid == inList.aminoAcid) ? this.totalDiff(inList) : this.totalCount() + inList.totalCount();
    }

    if(this.next == null || inList.next == null){

      if(this.next == null){
        //checks the possibility of a matching
        if (this.aminoAcid < inList.aminoAcid){
          return this.totalCount() + aminoCompareRes(inList);
        }if (this.aminoAcid > inList.aminoAcid){
          return inList.totalCount() + this.aminoAcidCompare(inList.next);
        }else {
          return this.totalDiff(inList) + aminoCompareRes(inList.next);
        }

      }else{
        //checks for a matching node
        if(inList.aminoAcid < this.aminoAcid ){
          return inList.totalCount() + aminoCompareRes(this);
        }if(inList.aminoAcid > this.aminoAcid){
          return this.totalCount() + this.next.aminoAcidCompare(inList);
        }else{
          return this.totalDiff(inList) + aminoCompareRes(this.next);
        }
      }
    }
    //if when amino matches
    if(this.aminoAcid == inList.aminoAcid){
      return this.totalDiff(inList) + this.next.aminoAcidCompare(inList.next);
    }

    if(this.aminoAcid < inList.aminoAcid){
      return this.totalCount() + this.next.aminoAcidCompare(inList);
    }

    return inList.totalCount() + this.aminoAcidCompare(inList.next);
  }

  /********************************************************************************************/
  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    //This is when they are in the last node
    if (this.next == null && inList.next == null){
      return (this.aminoAcid == inList.aminoAcid) ? this.codonDiff(inList) : this.totalCount() + inList.totalCount();
    }
    //last node refrence
    if(this.next == null || inList.next == null){
      //last node
      if(this.next == null){
        //checks matching
        if (this.aminoAcid < inList.aminoAcid){
          return this.totalCount() + aminoCompareRes(inList);
        }if (this.aminoAcid > inList.aminoAcid){
          return inList.totalCount() + this.aminoAcidCompare(inList.next);
        }else {
          return this.codonDiff(inList) + aminoCompareRes(inList.next);
        }

      }else {
        //checks matching
        if(inList.aminoAcid < this.aminoAcid ){
          return inList.totalCount() + aminoCompareRes(this);
        }if(inList.aminoAcid > this.aminoAcid){
          return this.totalCount() + this.next.aminoAcidCompare(inList);
        }else {
          return this.codonDiff(inList) + aminoCompareRes(this.next);
        }
      }
    }
    //ammino match
    if(this.aminoAcid == inList.aminoAcid){
      return this.codonDiff(inList) + this.next.codonCompare(inList.next);
    }
    //Cases when this has an aminoAcid that inList pointer does not
    if(this.aminoAcid < inList.aminoAcid){
      return this.totalCount() + this.next.codonCompare(inList);
    }
    //Cases when inList aminoAcid has an aminoacid that this does not
    return inList.totalCount() + this.codonCompare(inList.next);
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    // end of link list
    if(this.next == null){
      return new char[]{this.aminoAcid};
    }
    //Creates array
    char[] a = next.aminoAcidList();

    char[] r = new char[a.length+1];
    r[0] = aminoAcid;

    for (int i = 0; i < a.length; i++) {
      r[i+1] = a[i];
    }

    return r;
  }

  /********************************************************************************************/
  // rcusively returns
  public int[] aminoAcidCounts(){

    if(this.next == null){
      return new int[]{this.totalCount()};
    }

    int[] a = next.aminoAcidCounts();

    int[] r = new int[a.length+1];
    r[0] = this.totalCount();

    for (int i = 0; i < a.length; i++) {
      r[i+1] = a[i];
    }
    return r;

  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    //Base case
    if(this.next == null) return true;

    if(this.aminoAcid > this.next.aminoAcid){
      return false;
    }
    //Recursive call
    return this.next.isSorted();
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String sequence){

    if (sequence.length() < 3){
      return new AminoAcidLL();
    }
    // first codon
    String nextCodon = sequence.substring(0, 3);
    //Create the head
    AminoAcidLL head = null;

    while ((AminoAcidResources.getAminoAcidFromCodon(nextCodon) != '*') && (AminoAcidResources.getAminoAcidFromCodon(nextCodon) != (char)0)){

      if (head == null){
        head = new AminoAcidLL(nextCodon);
        sequence = sequence.substring(3);
        if(sequence.length() < 2) {
          nextCodon = sequence;
        }
        else {
          nextCodon = sequence.substring(0, 3);
        }

      }else {
        //Add a codon
        head.addCodon(nextCodon);
        sequence = sequence.substring(3);
        nextCodon = (sequence.length() == 3) ? sequence : sequence.substring(0,3);
      }
    }

    return head;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    if(inList == null){
      return new AminoAcidLL();
    }
//pointers

    AminoAcidLL beforeNode = inList;
    AminoAcidLL checkNode = beforeNode.next;
    AminoAcidLL nextNode = null;
    AminoAcidLL head = inList;
    AminoAcidLL insPos = null;
    //Loop till sorted
    while(checkNode != null){
      nextNode = checkNode.next;
      insPos = findInsPos(head, checkNode);
      if(insPos == beforeNode){
        beforeNode = checkNode;
      }
      else{
        beforeNode.next = checkNode.next;

        if(insPos == null){
          checkNode.next = head;
          head = checkNode;

        }
        else{
          checkNode.next = insPos.next;
          insPos.next = checkNode;
        }
      }
      checkNode = nextNode;
    }
    return head;
  }

  /********************************************************************************************/
  //method to see where to insert
  private static AminoAcidLL findInsPos(AminoAcidLL inList, AminoAcidLL aminoAcid){
    AminoAcidLL head = inList;
    char checkLetter = aminoAcid.aminoAcid;
    AminoAcidLL pos = null;
    AminoAcidLL checkPos = head;
    //Looks for after what node the node that is being checked should be introduced, where the aminoacid that is next has a greater value
    while(checkPos != null && checkLetter > checkPos.aminoAcid ){
      pos = checkPos;
      checkPos = checkPos.next;
    }
    return pos;
  }

  /********************************************************************************************/
  //this method helps check for the rest
  private static int aminoCompareRes(AminoAcidLL list){
    int count = 0;
    AminoAcidLL pointer = list;
    while(pointer != null){
      count += pointer.totalCount();
      pointer = pointer.next;
    }
    return count;
  }
}