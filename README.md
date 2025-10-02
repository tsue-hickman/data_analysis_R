# Overview

I created this R program to learn statistical computing and data visualization while working with real biological data. As someone interested in bioinformatics, I wanted to understand how researchers analyze gene expression patterns across different cancer types. This project helped me get comfortable with R's syntax, data structures, and visualization libraries.

The software reads normalized RNA-seq gene expression data and sample metadata, then performs variance-based analysis to identify genes with interesting expression patterns across five cancer types (LUAD, PRAD, BRCA, KIRC, COAD). It generates heatmaps to visualize these patterns and exports a list of high-variance genes for further investigation.

I wrote this to practice working with real-world datasets and to understand how bioinformaticians use R for exploratory data analysis. Plus, heatmaps look really cool and I wanted to make some.

Note: I had to completely delete and restart my original repository because my initial data.csv file was 196MB and exceeded GitHub's file size limits. Lesson learned about working with sample data first!

[Software Demo Video](http://youtube.link.goes.here)

# Development Environment

I used RStudio as my development environment, which made it easy to test code chunks and visualize outputs interactively. For version control, I used Git and GitHub Desktop, though I ran into some issues with large file management along the way.

The programming language is R (version 4.x). I used several libraries:

- **pheatmap** - for creating clustered heatmaps
- **RColorBrewer** - for color palettes in visualizations
- **genefilter** - for calculating row variances across gene expression data

The dataset format is CSV files with gene expression values (normalized log2-transformed data) and sample metadata indicating cancer types.

# Useful Websites

- [R Documentation](https://www.r-project.org/other-docs.html) - Official R documentation and manuals (specifically https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf)
- [pheatmap Documentation](https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf) - Helped me understand heatmap parameters
- [Stack Overflow - R Tag](https://stackoverflow.com/questions/tagged/r) - Answered specific syntax questions
- [CRAN Task Views](https://cran.r-project.org/web/views/) - Found relevant packages for bioinformatics
- [Quick R](https://www.statmethods.net/) - Good beginner-friendly R tutorials
- [Claude.Ai](https://claude.ai/) - Helped me plan the project and check over everything to ensure I hit all the benchmarks needed for the project
- [Grok.Ai](https://grok.com/) - Same purpose as Claude.Ai (helpful to have multiple perspectives)

# Future Work

- Add error handling for missing or malformed CSV files
- Implement additional visualization types (volcano plots, box plots by cancer type)
- Create functions to make the code more modular and reusable
- Add statistical tests to compare expression between cancer types
- Optimize performance for larger datasets
- Add command-line arguments to make file paths configurable
- Include data preprocessing steps for raw count data
