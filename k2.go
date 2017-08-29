package main

import (
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"fmt"
	"os"
	"strings"
	"math"
//	"reflect"
)

//score by probability
func probScorer(n, k, l int) (score float64) {
        //TODO
	return score
}

//score by sequence length
func lenScorer(n, k, l int) (score float64) {
	score = float64(n)/float64(l-k)
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
func maskSeq(sequ seq.Sequence, positions []int) (seqOut seq.Sequence) {
	firstpos := 0
	lastpos := 0
	seqOut = linear.NewSeq(sequ.Name(), []alphabet.Letter(""),alphabet.DNA)
	//fmt.Println(positions)
	for i := 0; i < len(positions)/2; i++ {
	//	fmt.Println(i)
	//	fmt.Println(seqOut)
		lastpos = positions[i*2]
	//	fmt.Println(lastpos)i
		seqOut = seqOut.AppendLetters(sequ.Slice().Slice(firstpos,lastpos))
	//	fmt.Println(seqOut)
		firstpos = lastpos
		lastpos = positions[i*2+1]
		seqOut = seqOut.AppendLetters(alphabet.Letter(strings.Repeat("N", lastpos - firstpos)))
        //        fmt.Println(seqOut)
	}
	seqOut = seqOut.Append(sequ.Slice().Slice(lastpos, sequ.Len()))
	return seqOut
}


func main() {
	//TODO command line input
	infile, err := os.Open("test.fa")
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer infile.Close()
	outfile, er := os.Create("testout_k2.fa")
	if er != nil {
                fmt.Println(err)
                os.Exit(1)
        }
	defer outfile.Close()
	//TODO allow choosing to filter or mask
	in := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
	out := fasta.NewWriter(outfile, 60)
	for {
		s, er := in.Read()

		if er != nil { //stop if EOF
			break

		} else { //process sequence
			//TODO variable window lengths
			winlen := 100
			k := 4
			//TODO variable k
			imap := kCounter(k,fmt.Sprintf("%v",s.Slice().Slice(0,winlen)))
			iscore := lenScorer(len(imap), k, winlen)
			fmt.Println(iscore)
			//TODO dynamic threshold
			threshold := 0.3
			filtering := false //toggle
			var filterpos []int = make([]int,0,10)
			if iscore < threshold {
				filtering = true
				filterpos = append(filterpos, 0)
			}
			for j := 0; j < s.Len() - winlen; j++ {
				if math.Mod(float64(j),1000000.0) == 0 {
					fmt.Println(j)
				}
				oldk := fmt.Sprintf("%v",s.Slice().Slice(j,j+k))
				newk := fmt.Sprintf("%v",s.Slice().Slice(j+winlen-k,j+winlen))
				if imap[oldk] == 1 {
					delete(imap, oldk)
				} else {
					imap[oldk]--
				}
				imap[newk]++
				iscore = lenScorer(len(imap), k, winlen)
				if filtering {
					if iscore > threshold {
						filterpos = append(filterpos, j + winlen)
						filtering = false
					}
				} else {
					if iscore < threshold {
						filterpos = append(filterpos, j)
						filtering = true
					}
				}
			}
			if filtering {
				filterpos = append(filterpos, s.Len())
			}
			if len(filterpos) == 0 {
				out.Write(s)
			} else {
				filtered := maskSeq(s, filterpos)
				out.Write(linear.NewSeq("",[]alphabet.Letter(fmt.Sprintf("%v",filtered.Slice(0,filtered.Len()))),alphabet.DNA))
			}
			//fmt.Println(imap)
			//fmt.Println(filtering)
			//fmt.Println(filterpos)
		}
	}
//	outfile.Close()
}

