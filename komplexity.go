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

func kQCounter(k int, seq string) (counts map[string]int) {
	counts = make(map[string]int)
	fmt.Println(len(seq))
	for i := 0; i < len(seq)-k; i++ {
		fmt.Println(seq[i:i+k])
		counts[seq[i:i+k]]++
	}
	fmt.Println(counts)
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
func maskQSeq(seqi string, sequ seq.Sequence, positions []int) (seqOut *linear.QSeq) {
	beghq := 0
	endhq := 0
	seqRaw := "" //sequ.Slice().Slice(0,0)

	//var qlett *alphabet.QLetters

	//iteratively mask by replacing sequence with N
	for i := 0; i < len(positions)/2; i++ {
		endhq = positions[i*2]
		seqRaw += seqi[beghq:endhq]
		beghq = positions[i*2+1]
		seqRaw += strings.Repeat("N", beghq - endhq)
		//seqOut = seqOut.Append(linear.NewQSeq("",[]alphabet.QLetter(strings.Repeat("N", beghq - endhq)),alphabet.DNA,alphabet.Sanger).Slice())
	}
	seqRaw += seqi[beghq:len(seqi)]
	//fmt.Println(s.Slice().(alphabet.QLetters)))
	//seqOut = seqOut.Append(sequ.Slice().Slice(beghq, sequ.Len()))
	seqOut = linear.NewQSeq(sequ.Name(), nil, alphabet.DNA, alphabet.Sanger)
	seqRawLetters := []alphabet.Letter(seqRaw)
	for j := range len(seqi) {
		seqOut.AppendQLetters(alphabet.QLetter{L: seqRawLetters[j], Q: 0})
	}
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
		//fmt.Printf("%T\n", s)
		//fmt.Printf("%T\n", in)
		//fmt.Printf("%-s\n", s)
		//fmt.Printf("%-s\n", in)
		qq := fmt.Sprint(s.Slice())
		fmt.Println(qq[1:len(qq)-1])
		// outfq.Write(s)
		//create map
		// var i_slice string
		var imap map[string]int
		fqs := ""
		if format == "fa" {
			a_slice := s.Slice().Slice(0,winlen)
			fmt.Printf("%T\n", a_slice)
			imap = kCounter(k,fmt.Sprintf("%v",a_slice))
		} else {
			b_slice := s.Slice().Slice(0,winlen)
			fmt.Printf("%T\n", b_slice)
			fmt.Printf("%-s\n", b_slice)
			fmt.Println(strings.Replace(fmt.Sprint(b_slice)[1:(winlen*2)], " ", "", -1))
			fqs = strings.Replace(fmt.Sprint(s.Slice().Slice(0,s.Len()))[1:(s.Len()*2)], " ", "", -1)
			fmt.Println(fqs)
			imap = kCounter(k,strings.Replace(fmt.Sprint(b_slice)[1:(winlen*2)], " ", "", -1))
		}


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

			var oldk string
			var newk string
			//grab first and last kmers
			if format == "fa" {
				oldk = fmt.Sprintf("%v",s.Slice().Slice(j,j+k))
				newk = fmt.Sprintf("%v",s.Slice().Slice(j+winlen-k,j+winlen))
			} else {
				oldk = fqs[j:j+k]
				newk = fqs[j+winlen-k:j+winlen]
			}
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
		fmt.Println(filterpos)
		if len(filterpos) == 0 {
			if format == "fa" {
				outfa.Write(s)
			} else {
				outfq.Write(s)
			}
		} else {
			if format == "fa" {
				fafiltered := maskSeq(s, filterpos)
				//fmt.Println(fafiltered)
				outfa.Write(linear.NewSeq(s.Name(),[]alphabet.Letter(fmt.Sprintf("%v",fafiltered.Slice(0,fafiltered.Len()))),alphabet.DNA))
			} else {
				fqfiltered := maskQSeq(fqs, s, filterpos)
				fmt.Println("up")
				fmt.Println(fqfiltered)
				fmt.Println("down")
				outfq.Write(fqfiltered)
				//outfq.Write(linear.NewQSeq(s.name(), []alphabet.QLetter))
			}
		}
	}
}
