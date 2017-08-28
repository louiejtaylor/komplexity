package main

import (
//	"github.com/biogo/biogo/complexity"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
//	"github.com/biogo/biogo/seq"
	"fmt"
	"os"
//	"strings"

//	"strconv"

	"reflect"
)

//score by probability
func probScorer(n, k, l int) (score int) {
        return score
}

//score by sequence length
func lenScorer(n, k, l int) (score int) {
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
	return seqOut
}


func main() {
	//in = fasta.NewReader("test.fa", linear.NewSeq("", nil, alphabet.DNA))
	//var sequ = linear.NewSeq("test_seq",[]alphabet.Letter("AGTCGATCAGCTAGCCGACTAGCAC"), alphabet.DNA)
	//fmt.Println(reflect.TypeOf(sequ.Features))

	//var in = *fasta.Reader

	infile, err := os.Open("test.fa")
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	in := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
	for {
		z, er := in.Read()
		if er != nil {
			fmt.Println(err)
			break//panic("!!")
		} else {
			y := fmt.Sprintf("%v",z.Slice().Slice(1,10))// strings.Join(z.Slice(),"")//.(alphabet.Letters)
			//fmt.Println(kCounter(4,z))
			fmt.Println(reflect.TypeOf(y))

			fmt.Println(y)
			fmt.Println(kCounter(4,y))
//			x := reflect.TypeOf(y).Elem()
//			InspectStruct(y)
//			for j := 0; j<x.NumField(); j++ {
//				fmt.Println(x.Field(j).Name)
//			}
//			fmt.Println(reflect.ValueOf(y))
		}
	}

	//var z, e = complexity.Z(sequ,sequ.Start(),sequ.End())
	//fmt.Println(z)
	//fmt.Println(e)

	//testing kmer counter
	fmt.Println(kCounter(4,"AAAAACCCC"))
	//panic("AAAH")
/*
	//TESTING WHETHER A ZERO MAP VALUE IS INCLUDED IN len
	//spoiler: it isn't
	fmt.Println("--tests--")
	tm := make(map[string]int)
	fmt.Println(len(tm))
	tm["test"] = tm["test"] + 1
	fmt.Println(len(tm))
	tm["test"] = tm["test"] - 1
	fmt.Println(len(tm))
	delete(tm,"test")
	fmt.Println(len(tm))

*/
}

