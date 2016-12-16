# BioPatML FSharp

This project is the work of a computer science student at Queensland University of Technology, under the supervision of Associate Professor James Hogan. The aim of the project was an exploration of using FSharp with the .NET Bio Bioinformatics library, and the eventual development of a Type Provider for GenBank data.

The goal here was fleshed out to involve implementing a workflow for identifying certain traits within a genome which would eventually be facilitated by the type provider. It was also decided that it would be worthwhile to take advantage of the powerful pattern matching capabilities of the F# language to facilitate pattern identification in the context of gene sequences. Particular attention was brought to a project developed by Stefan Maetschke by the name of “BioPatML” and how functionality similar to that of BioPatML would be of significant benefit to the project. As such, the primary aim of my part the project became to implement a version of the existing BioPatML project in the F# language

## BioPatML
BioPatML is a pattern matching language that allows users to define complex patterns that exist within gene sequences. It aims to provide enhanced pattern description capabilities that cannot be achieved accu- rately with regular expressions or position weight matrices. The language is XML based, whereby the user can describe patterns as individual XML components as well as more structured patterns which contain a collection of these smaller patterns in some structured context.

### Patterns
BioPatML supports several different pattern variants. These can be broadly categorised under each of the following:

* Region Patterns
    Describe a region within a DNA sequence without being specific about its contents.
* Motif Patterns
    Describe the specific contents of a pattern (such as a certain set of nucleotides, e.g. “ATTG”).
* Structured Patterns
    Describe patterns which are aggregations of several lower-level patterns (such as Motif patterns).
    
A full outline of the range of supported pattern types can be found in the Wiki, accompanied by examples.

# Known Issues and Future Directions
The BioPatML F# library contains versions of most, but not all patterns existing within the original version of BioPatML. The specifics of this support can be seen below:

| Pattern | Notes | Supported | 
| --- | --- | --- |
| Any | Full Support | Y |
| Gap | Full Support | Y | 
| Position-Weight Matrix | - | N |
| Anchor | - | N | 
| Motif | Full Support | Y | 
| Regex | Full Support | Y | 
| Prosite | Support with Limitations | Y |
| Repeat | Support with Limitations | Y |
| Set | Support with Limitations | Y | 
| Series | Support with Limitations | Y | 
| Profile | - | N | 

## Future Work 

Some future work that can be done to improve upon the current version of the BioPatML F# library includes:

* Implementing support for pattern types not currently supported.
    Position Weight Matrix – Anchor
    Profile
* Re-evaluating the type architecture and the multitude of issues arising from circular dependency-related issues.
* Implementing further features for supported types that currently have limitations.
    Prosite
    
    --Revised matching model
    
    --Support for repeat ranges
      
    Series and Repeat
    
    --Support for MatchablePat types other than Motif
    
    Set
    
    --Implementing additional match process that performs matches for every component and spec- ifies that with the highest score
    
    --Potential support for other Set patterns
      
* More rigorous and automated testing to ensure that no bugs have slipped under the radar.
