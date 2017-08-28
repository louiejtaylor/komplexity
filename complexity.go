package main

import (
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"fmt"
	"os"
)

//score by probability
func probScorer(n, k, l int) (score int) {
        //TODO
	return score
}

//score by sequence length
func lenScorer(n, k, l int) (score int) {
	//TODO
        return score
}

//return kmers and frequencies in a given sequence
func kCounter(k int, seq string) (counts map[string]int) {
	counts = make(map[string]int)
	for i := 0; i < len(seq)-k; i++ {
		counts[seq[i:i+k]]++
	}
        return counts
}
//mask a sequence given a list of positions [start, end, start, end...]
func maskSeq(seq string, positions []int) (seqOut string) {
	//TODO
	return seqOut
}


func main() {
	//TODO command line input
	infile, err := os.Open("test.fa")
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	//TODO allow choosing to filter or mask
	in := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
	for {
		z, er := in.Read()
		if er != nil { //stop if reading error
			fmt.Println(err)
			break
		} else { //process sequence
			//TODO variable window lengths
			winlen := 100
			for j := 0; j < z.Len() - winlen+1; j++ {
				y := fmt.Sprintf("%v",z.Slice().Slice(j,j+winlen+1))

				fmt.Println(y)
				fmt.Println(kCounter(4,y))
			}
		}
	}
}

