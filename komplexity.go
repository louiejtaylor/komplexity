package main

import (
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq"

	"fmt"
	"os"
	"strings"
	"flag"
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

	//iteratively mask by replacing sequence with N
	for i := 0; i < len(positions)/2; i++ {
		endhq = positions[i*2]
		seqOut = seqOut.Append(sequ.Slice().Slice(beghq,endhq))
		beghq = positions[i*2+1]
		seqOut = seqOut.Append(linear.NewSeq("",[]alphabet.Letter(strings.Repeat("N", beghq - endhq)),alphabet.DNA).Slice())
	}
	seqOut = seqOut.Append(sequ.Slice().Slice(beghq, sequ.Len()))
	return seqOut
}

//mask a sequence with quality scores // this should probably be merged
func maskQSeq(sequ seq.Sequence, positions []int) (seqOut alphabet.Slice) {
	beghq := 0
	endhq := 0
	seqOut = sequ.Slice().Slice(0,0)
	//var qlett *alphabet.QLetters

	//iteratively mask by replacing sequence with N
	for i := 0; i < len(positions)/2; i++ {
		endhq = positions[i*2]
		seqOut = seqOut.Append(sequ.Slice().Slice(beghq,endhq))
		beghq = positions[i*2+1]

		//seqOut = seqOut.Append(linear.NewQSeq("",[]alphabet.QLetter(strings.Repeat("N", beghq - endhq)),alphabet.DNA,alphabet.Sanger).Slice())
	}
	seqOut = seqOut.Append(sequ.Slice().Slice(beghq, sequ.Len()))
	return seqOut
}

func main() {
	//grab input from command line
	var fa string
	flag.StringVar(&fa, "fa", "", "Input .fasta file")

	var fq string
	flag.StringVar(&fq, "fq", "", "Input .fastq file")

	var k int
	flag.IntVar(&k, "k", 4, "k-mer size")

	var winlen int
	flag.IntVar(&winlen,"win", 100, "window length")

	var fileout string
	flag.StringVar(&fileout, "out", "", "Output filename (default: appends '_filtered' to input fname")

	flag.Parse()

	format := ""
	filein := ""

	if fa == "" && fq == "" {
		fmt.Println("Must provide either -fq or -fa file")
		fmt.Println(flag.Usage)
	} else if fa != "" && fq != "" {
		fmt.Println("Must provide only one of -fq or -fa")
		fmt.Println(flag.Usage)
	} else {
		if fa == "" {
			format = "fq"
			filein = fq
		} else {
			format = "fa"
			filein = fa
		}
	}

	//set up files
	if fileout == "" {
		fileout = filein+"_filtered"
	}

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

	var in *seqio.Scanner
	var outfa *fasta.Writer
	var outfq *fastq.Writer
	//setup reader and writer
	//in := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
	if format == "fa" {
		in = seqio.NewScanner(fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA)))
		outfa = fasta.NewWriter(outfile, 60)
	} else {
		in = seqio.NewScanner(fastq.NewReader(infile, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger)))
		outfq = fastq.NewWriter(outfile)

		//fmt.Println("fastq handling not yet implemented")
		//os.Exit(1)
	}
	//read infile
	for in.Next() {
		s := in.Seq()//.(*linear.QSeq)//.(*linear.Seq) //does QSeq work here?
		fmt.Println(s.Name())
		fmt.Println(s.Slice().Slice(0,winlen))
		// outfq.Write(s)
		//create map
		imap := kCounter(k,fmt.Sprintf("%v",s.Slice().Slice(0,winlen)))
		iscore := lenScorer(len(imap), k, winlen)

		threshold := 0.55 //make configurable
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
		// fmt.Println(filterpos)
		//after processing, add final position
		if filtering {
			filterpos = append(filterpos, s.Len())
		}
		if len(filterpos) == 0 {
			if format == "fa" {
				outfa.Write(s)
			} else {
				outfq.Write(s)
			}
		} else {
			if format == "fa" {
				fafiltered := maskSeq(s, filterpos)
				fmt.Println(fafiltered)
				outfa.Write(linear.NewSeq(s.Name(),[]alphabet.Letter(fmt.Sprintf("%v",fafiltered.Slice(0,fafiltered.Len()))),alphabet.DNA))
			} else {
				fqfiltered := maskQSeq(s, filterpos)
				fmt.Println("up")
				fmt.Println(fqfiltered)
				fmt.Println("down")
				outfq.Write(s)
				//outfq.Write(linear.NewQSeq(s.name(), []alphabet.QLetter))
			}
		}
	}
}
