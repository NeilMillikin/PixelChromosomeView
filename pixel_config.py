
# There are only TWO MANDATORY settings:
DATA_FILE_DIRECTORY = 'demo_raw_dna'
siblings_to_render = ['JULIE', 'ALLISON', 'COLLETTE']

"""
NOTE: raw DNA file names must:
1) contain name of source vendor (i.e. Ancestry, 23andMe, MyHeritage, etc
2) contain person's name exactly as listed in 'siblings_to_render' below
3) contain the word 'raw'
4) lastly, raw data musst must be in  .txt files

example: 23andMe_JULIE_raw_dna.txt

copy all relevant raw dna files into the directory your specify
"""
# optionally you can compare the siblings with one additional relative
extra_match = ''

# removes a lot of 'noise' SNPs that don't contribute to the analysis
FILTER_COMPLETELY_MATCHED_SEGMENTS = True

##### optional settings to change appearance of rendered chromosome pairs  ###
FULLY_IDENTICAL_SNP_COLOR = 'limegreen'
NO_MATCH_SNP_COLOR = 'crimson'
HALF_IDENTICAL_SNP_COLOR = 'yellow'
BACKGROUND_COLOR = 'lightsteelblue'
CHROMOSOME_BASE_COLOR = 'white'
WIDTH_OF_SNP_LINE = 1
HEIGHT_OF_CHROMOSOME_IMAGE = 100
SPACE_BETWEEN_MATCHES = 50
CHROM_PAGE_LEFT_BORDER = 450
CHROM_PAGE_TEXT_BORDER = 50
CHROM_TITLE_TEXT_FONT_SIZE = 60
CHROM_MATCH_TEXT_FONT_SIZE = 40
CHROM_PAGE_RIGHT_BORDER = 50
CHROM_PAGE_TOP_BORDER = 150
CHROM_PAGE_TITLE_SPACE = 200
CHROM_PAGE_BOTTOM_BORDER = 150
MINIMUM_PAGE_WIDTH = 1750
MILESTONE_SPACING = 50
TICKER_SPACING = 10
MILESTONE_VERTICAL_POSITION = -25
TICKER_VERTICAL_POSITION = -15
TICKMARK_FONT = "arial.tff"
TICKMARK_FONT_SIZE = 10

