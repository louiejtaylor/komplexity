package main

import (
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"

	"fmt"
	"os"
	"strings"
//	"math"
	"flag"

//	"reflect"
)

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
func maskSeq(sequ seq.Sequence, positions []int) (seqOut alphabet.Slice) {
	beghq := 0
	endhq := 0
	seqOut = sequ.Slice().Slice(0,0)
	for i := 0; i < len(positions)/2; i++ {
		endhq = positions[i*2]
		seqOut = seqOut.Append(sequ.Slice().Slice(beghq,endhq))
		beghq = positions[i*2+1]
		seqOut = seqOut.Append(linear.NewSeq("",[]alphabet.Letter(strings.Repeat("N", beghq - endhq)),alphabet.DNA).Slice())
	}
	seqOut = seqOut.Append(sequ.Slice().Slice(beghq, sequ.Len()))
	return seqOut
}


func main() {
	//grab input from command line
	var filein string
	flag.StringVar(&filein, "in", "test.fa", "Input filename")

	//k := flag.Int("k", 4, "k-mer size")
	var k int
	flag.IntVar(&k, "k", 4, "k-mer size")

	//winlen := flag.Int("win", 100, "window length")
	var winlen int
	flag.IntVar(&winlen,"win", 100, "window length")

	var fileout string
	flag.StringVar(&fileout, "out", filein+"_filtered", "Output filename")

	flag.Parse()

	//set up files
	infile, err := os.Open(filein)
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	defer infile.Close()
	outfile, er := os.Create(fileout)
	if er != nil {
                fmt.Println(err)
                os.Exit(1)
        }
	defer outfile.Close()

	//setup reader and writer
	in := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
	out := fasta.NewWriter(outfile, 60)

	//read infile
	for {
		s, er := in.Read()

		if er != nil { //stop if EOF
			break

		} else { //process sequence

			imap := kCounter(k,fmt.Sprintf("%v",s.Slice().Slice(0,winlen)))
			iscore := lenScorer(len(imap), k, winlen)
			fmt.Println(iscore)

			threshold := 0.55
			filtering := false //toggle

			var filterpos []int = make([]int,0,10)

			//check if initial score below threshold
			if iscore < threshold {
				filtering = true
				filterpos = append(filterpos, 0)
			}

			for j := 0; j < s.Len() - winlen; j++ {

				//grab first and last kmers
				oldk := fmt.Sprintf("%v",s.Slice().Slice(j,j+k))
				newk := fmt.Sprintf("%v",s.Slice().Slice(j+winlen-k,j+winlen))

				//get rid of old kmer
				if imap[oldk] == 1 {
					delete(imap, oldk)
				} else {
					imap[oldk]--
				}

				//add new kmer
				imap[newk]++
				iscore = lenScorer(len(imap), k, winlen)

				//add "inflection point" if necessary
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

			//after processing, add final position
			if filtering {
				filterpos = append(filterpos, s.Len())
			}

			if len(filterpos) == 0 {
				out.Write(s)
			} else {
				filtered := maskSeq(s, filterpos)
				out.Write(linear.NewSeq("ab",[]alphabet.Letter(fmt.Sprintf("%v",filtered.Slice(0,filtered.Len()))),alphabet.DNA))
			}
		}
	}
}

