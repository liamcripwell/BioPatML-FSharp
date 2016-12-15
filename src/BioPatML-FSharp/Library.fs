namespace Bio.FSharp

//          +-------------------------------+          //
//          |  BioPatML: F# Implementation  |          //
//----------|                               |----------//
//          |     Author: Liam Cripwell     |          //
//          +-------------------------------+          //

module BioPatML =
    open Bio.Core
    open FSharp.Data
    open System
    open System.IO
    open System.Collections
    open Bio.Core.Extensions
    open Bio.Extensions
    open System.Text.RegularExpressions

    //          +---------------------+          //
    //----------|  Utility Functions  |----------//
    //          +---------------------+          //

    let (|RegexPattern|_|) pattern input =
       let m = Regex.Match(input,pattern) 
       if (m.Success) then Some m.Groups.[1].Value else None 

    type Alphabet =
    | DNA
    | RNA
    | PROTEIN

    // checks whether a char is a valid DNA alphabet member
    let isDNA (item : char) =
        Bio.DnaAlphabet.Instance.GetValidSymbols().Contains(Convert.ToByte(item))

    let isRNA (item: char) =
        Bio.RnaAlphabet.Instance.GetValidSymbols().Contains(Convert.ToByte(item))

    let isProtein (item: char) =
        Bio.ProteinAlphabet.Instance.GetValidSymbols().Contains(Convert.ToByte(item))

    // determines whether a given value is a wildcard value
    let isWildcard chunk =
        match chunk.ToString().ToLower() with
        | x when x = "x" -> true
        | x when x = "n" -> true
        | _ -> false

    // checks whether a list of chars are all valid DNA alphabet members
    let checkAlphabet alphabet slice =
        let validator =
                match alphabet with
                | DNA     -> isDNA
                | RNA     -> isRNA
                | PROTEIN -> isProtein

        let wildOrValid nuc =
            match nuc with
            | x when validator x  -> true
            | x when isWildcard x -> true
            | _                   -> false

        (List.filter wildOrValid slice).Length = slice.Length

                     
    //          +---------------------+          //
    //----------|  Abstract Patterns  |----------//
    //          +---------------------+          //

    [<AbstractClass>]
    type BioPat = 
        new () = { }

    [<AbstractClass>]
    type RegionPat(minLength, maxLength) =
        inherit BioPat()
        member val Min = minLength with get
        member val Max = maxLength with get

    [<AbstractClass>]
    type MatchablePat() =
        inherit BioPat()
        abstract member Match : Bio.ISequence -> bool

    [<AbstractClass>]
    type MotifPat(pattern:String) =
        inherit MatchablePat()
    
    [<AbstractClass>]
    type StructuredPat(content:List<BioPat>) =
        inherit MatchablePat()


    //          +-------------------+          //
    //----------|  Region Patterns  |----------//
    //          +-------------------+          //

    /// A pattern that matches any sequence within minLength and maxLength
    type Any(minLength, maxLength) =
        inherit RegionPat(minLength, maxLength)

        member this.Match(input:Bio.ISequence) =
            let inputString = input.ConvertToString()

            let rec examine input acc i =
                match input with
                | x::xs -> if List.length input >= i && i <= maxLength then 
                               examine input (List.append acc [(List.take i input)]) (i+1)
                           else 
                               examine xs acc minLength
                | _     -> acc

            examine (List.ofSeq inputString) List.empty minLength

    /// Represents a gap within existing between two patterns
    type Gap(minLength, maxLength) =
        inherit RegionPat(minLength, maxLength)


    //          +------------------+          //
    //----------|  Motif Patterns  |----------//
    //          +------------------+          //

    /// A pattern that finds matches of a short sequence above a similarity thershold
    type Motif(pattern:String, ?threshold, ?alphabet) =
        inherit MotifPat(pattern)

        let threshold      = defaultArg threshold 1.0
        let alphabet       = defaultArg alphabet DNA
        let motif          = pattern.ToLower() |> Seq.toList
        let weightFraction = 1.0 / float motif.Length

        // will fail construction if motif is not a valid DNA sequence
        do
            if not (checkAlphabet alphabet motif) then 
                failwith "Supplied motif is not valid DNA sequence"

        override this.Match(input:Bio.ISequence) =
            let inputString = input.ConvertToString().ToLower()
            let motifString = motif |> Array.ofList |> System.String

            let rec examine motifSlice inputSlice =
                match motifSlice with
                | x::xs when xs <> []
                        -> let isMatch = x |> (fun symbol ->
                               match symbol with
                               | x when isWildcard x   -> true
                               | x when x = 
                                   Seq.head inputSlice -> true
                               | _                     -> false);
                           match isMatch with
                           | true ->
                               // handle case where inputSlice has matched completely but for its 
                               // length being shorter than the motif
                               match (inputSlice |> Seq.tail |> List.ofSeq) with
                               | []  -> float (Seq.length motif - Seq.length motifSlice) * weightFraction
                               | _   -> examine xs (Seq.tail inputSlice)
                           | false ->
                               float (Seq.length motif - Seq.length motifSlice) * weightFraction
                | _     -> 1.0
        
            match this.Threshold with
                | 1.0 when (inputString.Length = motifString.Length) -> inputString = motifString
                | _   -> examine motif inputString >= this.Threshold
            |> (fun x -> this.Threshold <- threshold; x) // reset threshold back to original value

        member val Threshold = threshold with get, set
        member val Motif     = motif with get

    /// Pattern that matches input against a regular expression
    type BioRegex(pattern:String) =
        inherit MotifPat(pattern)

        override this.Match(input) =
            let inputString = input.ConvertToString()
            let m = Regex.Match(inputString.ToLower(), pattern.ToLower())
            m.Success

    /// Pattern that matches input agaisnt a structured repeating sequence of motifs
    type Repeat(pattern:String, minGap, maxGap, ?threshold, ?repeatCount) =
        inherit MotifPat(pattern)

        let threshold = defaultArg threshold 1.0
        let repeatCount = defaultArg repeatCount 1
        let gap = Gap(minGap, maxGap)
        let motif = Motif(pattern, threshold)

        let rec constructComponents (acc: List<BioPat>) i =
            match i with
            | x when x < repeatCount -> 
                constructComponents (List.append acc [gap; motif]) (i+1)
            | _ -> acc

        let components = constructComponents [motif] 0

        let gaps   = components |> List.filter (fun x -> (x.GetType() = typeof<Gap>)) 
                                |> List.map (fun x -> x :?> Gap)
        let motifs = components |> List.filter (fun x -> (x.GetType() = typeof<Motif>)) 
                                |> List.map (fun x -> x :?> Motif)

        override this.Match(input:Bio.ISequence) = 
            let splitPoint motRef gapSize = 
                motifs.Item motRef |> (fun x -> x.Motif.Length)
                                   |> (+) gapSize

            let splitDNA (dna: List<'a>) motRef gapSize =
                let point = splitPoint motRef gapSize
                match point with
                | x when x >= dna.Length -> None
                | _ -> Some <| List.splitAt point dna

            // finds all possible gaps at a point in dna seq given last motif index and remaining dna slice
            let getGaps mot slice =
                slice |> List.ofSeq
                      |> (fun x -> List.map (fun gapSize -> splitDNA x mot gapSize) [(gaps.Item mot).Min .. (gaps.Item mot).Max])
                      |> List.filter Option.isSome
                      |> List.map (fun x -> snd <| Option.get x)

            // determines if the input matches the Series
            let isMatch input = 
                let rec compute mot (slice: String) = 
                    if (motifs.Item mot).Match(new Bio.Sequence(Bio.Alphabets.DNA, slice)) then 
                        if mot < gaps.Length then
                            slice |> getGaps mot
                                  |> List.map (fun x -> x |> List.toArray |> (fun s -> System.String s))
                                  |> List.map (compute (mot+1))
                                  |> List.contains true
                        else
                            true // no more gaps and last motif was matched
                    else 
                       false

                compute 0 input // begin computation with first motif

            // return boolean result of match
            input.ConvertToString() |> isMatch

    /// Prosite Pattern
    type Prosite(pattern:String) = 
        inherit MotifPat(pattern)

        let parseRegex token = 
            token |> Some

        let parseExclusion token =
            match Regex.Match(token, "\{[a-zA-Z]+\}").Success with
            | true  -> let inner = token.Trim([| '{'; '}' |]);
                       match inner |> List.ofSeq 
                                   |> (checkAlphabet DNA) with
                       | true  -> String.concat "" ["[^"; inner; "]"]
                                   |> Some
                       | false -> None
            | false -> None

        let parseRepeat token =
            match Regex.Match(token, "[a-zA-Z]\([0-9]+,*[0-9]+\)").Success with
            | true  -> let content = (token.Split [|'('|]).[0];
                       match content |> List.ofSeq 
                                     |> (checkAlphabet DNA) with
                       | true  -> Some <| token.Replace("(", "{").Replace(")", "}")
                       | false -> None
            | false -> None

        let parseWildcard =
            "." |> Some

        let parseGeneral token =
            token |> Seq.toList
                  |> (checkAlphabet DNA)
                  |> (fun res -> match res with
                                 | true -> token |> Some
                                 | false -> None)

        let parseToken (token:String): Option<string> =
            match token with
            | x when x.Contains("[") -> parseRegex x
            | x when x.Contains("{") -> parseExclusion x
            | x when x.Contains("(") -> parseRepeat x
            | x when isWildcard x    -> parseWildcard
            | x                      -> parseGeneral x

        let tokens = pattern.Split[|'-'|]

        let pattern =
            tokens |> List.ofArray
                   |> List.map parseToken
                   |> (fun res -> match res with
                                  | x when List.contains None x -> 
                                               failwith "Invalid Prosite pattern"
                                  | x -> x |> List.map Option.get)
                   |> List.fold (+) ""

        override this.Match(input) =
            BioRegex(pattern).Match(input)


    //          +-----------------------+          //
    //----------|  Structured Patterns  |----------//
    //          +-----------------------+          //

    /// A Series of patterns which aggregate to form a new, more nuanced pattern
    type Series(content : List<BioPat>) =
        inherit StructuredPat(content)

        // seperate Motifs and Gaps into their own lists
        let gaps   = content |> List.filter (fun x -> (x.GetType() = typeof<Gap>)) 
                             |> List.map (fun x -> x :?> Gap)
        let motifs = content |> List.filter (fun x -> (x.GetType() = typeof<Motif>)) 
                             |> List.map (fun x -> x :?> Motif)
 
        override this.Match(input:Bio.ISequence) = 
            let splitPoint motRef gapSize = 
                motifs.Item motRef |> (fun x -> x.Motif.Length)
                                   |> (+) gapSize

            let splitDNA (dna: List<'a>) motRef gapSize =
                let point = splitPoint motRef gapSize
                match point with
                | x when x >= dna.Length -> None
                | _ -> Some <| List.splitAt point dna

            // finds all possible gaps at a point in dna seq given last motif index and remaining dna slice
            let getGaps mot slice =
                slice |> List.ofSeq
                      |> (fun x -> List.map (fun gapSize -> splitDNA x mot gapSize) [(gaps.Item mot).Min .. (gaps.Item mot).Max])
                      |> List.filter Option.isSome
                      |> List.map (fun x -> snd <| Option.get x)

            // determines if the input matches the Series
            let isMatch input = 
                let rec compute mot (slice: String) = 
                    if (motifs.Item mot).Match(new Bio.Sequence(Bio.Alphabets.DNA, slice)) then 
                        if mot < gaps.Length then
                            slice |> getGaps mot
                                  |> List.map (fun x -> x |> List.toArray |> (fun s -> System.String s))
                                  |> List.map (compute (mot+1))
                                  |> List.contains true
                        else
                            true // no more gaps and last motif was matched
                    else 
                       false

                compute 0 input // begin computation with first motif

            // return boolean result of match
            input.ConvertToString() |> isMatch

    type SettableUnion =
            | M of Motif
            | X of BioRegex
            | R of Repeat
            | P of Prosite
            | S of Series

    /// Can be used to house a collection of patterns which can be matched against
    type Set(content:List<BioPat>, ?threshold) =
        inherit StructuredPat(content)

        let threshold = defaultArg threshold 1.0

        // temporarily overrides a Motif's match threshold
        let updatePush (motif: Motif) = (motif.Threshold <- threshold); motif

        let getPatternType (pattern: BioPat) =
            match pattern with
            | x when x.GetType() = typeof<Motif>    -> (x :?> Motif)    |> updatePush 
                                                                        |> M |> Some
            | x when x.GetType() = typeof<BioRegex> -> (x :?> BioRegex) |> X |> Some
            | x when x.GetType() = typeof<Repeat>   -> (x :?> Repeat)   |> R |> Some
            | x when x.GetType() = typeof<Prosite>  -> (x :?> Prosite)  |> P |> Some
            | x when x.GetType() = typeof<Series>   -> (x :?> Series)   |> S |> Some
            | _                                     -> None

        // extracts Match function from pattern and then performs it on input
        let doMatch input pattern =
            match pattern |> Option.isSome with
            | false -> failwith "Supplied pattern is invalid"
            | true  -> pattern 
                       |> Option.get
                       |> (fun x -> 
                               match x with
                               | M pat -> pat.Match
                               | X pat -> pat.Match
                               | R pat -> pat.Match
                               | P pat -> pat.Match
                               | S pat -> pat.Match)
                       |> (fun matcher -> matcher(input))

        // gracefully handles binding of function to (input, pattern) tuple
        let (>>=) (input, pattern) f =
            let castPat = getPatternType pattern
            f input castPat

        override this.Match(input) =
            let inputString = input.ConvertToString()

            content |> List.tryPick (fun x -> match (input, x) >>= doMatch with
                                              | true  -> Some true
                                              | false -> None)
                    |> (fun x -> match x with
                                 | Some y -> Option.get x
                                 | None   -> false)


    //          +-------------------------+          //
    //----------|  User-Facing Functions  |----------//
    //          +-------------------------+          //

    type MatchableUnion =
        | M of Motif
        | X of BioRegex
        | R of Repeat
        | P of Prosite
        | S of Set
        | Z of Series

    let getPatternType (pattern: MatchablePat) =
        match pattern with
        | x when x.GetType() = typeof<Motif>    -> Some <| (M <| (x :?> Motif))
        | x when x.GetType() = typeof<BioRegex> -> Some <| (X <| (x :?> BioRegex))
        | x when x.GetType() = typeof<Repeat>   -> Some <| (R <| (x :?> Repeat))
        | x when x.GetType() = typeof<Prosite>  -> Some <| (P <| (x :?> Prosite))
        | x when x.GetType() = typeof<Set>      -> Some <| (S <| (x :?> Set))
        | x when x.GetType() = typeof<Series>   -> Some <| (Z <| (x :?> Series))
        | _                                     -> None

    let doMatch input pattern =
        match pattern |> Option.isSome with
        | false -> failwith "Supplied pattern is invalid"
        | true  -> pattern 
                   |> Option.get
                   |> (fun x -> 
                           match x with
                           | M pat -> pat.Match
                           | X pat -> pat.Match
                           | R pat -> pat.Match
                           | P pat -> pat.Match
                           | S pat -> pat.Match
                           | Z pat -> pat.Match)
                   |> (fun matcher -> matcher(input))

    let (>>=) (input, pattern) f =
        let castPat = getPatternType pattern
        f input castPat

    /// Locates first occurence of pattern within a dna sequence
    let locate (input: Bio.ISequence) pattern = 
        let inputString = input.ConvertToString()

        let rec compute slice =
            match slice with
            | x::xs -> let bio = (Bio.Sequence(Bio.Alphabets.DNA, slice 
                                               |> Array.ofList 
                                               |> System.String))
                       match (bio, pattern) >>= doMatch with 
                       | false -> compute slice.Tail
                       | true  -> Some <| inputString.Length - slice.Length
            | _     -> None

        inputString |> List.ofSeq |> compute

    /// Checks whether a pattern exists within a dna sequence
    let exists input pattern =
        locate input pattern |> Option.isSome