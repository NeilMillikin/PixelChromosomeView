#!/usr/bin/env python3
"""
    Copyright 2021 by Neil Millikin
    contact: neil.millikin@gmail.com
    license: GPLv3

    This program graphically displays DNA match data for three or more siblings

    The comparison of three or more siblings' DNA data has been shoen to allow
    derivation of the maternal and paternal grandparents for each sibling
    chromosome pair.

    The process for analyzing similar data sets has been termed "Visual Phasing"

    To use this program,  original raw DNA must be downloaded from
    an original DNA testing vendor, such as AncestryDNA, MyHeritage, 23andMe
    or similar sources.  Once inserted into this program,
    the data is parsed, analyzed and graphically represented.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

__author__ = "Neil Millikin"
__maintainer__ = "Neil Millikin"
__contact__ = "neil.millikin@gmail.com"
__email__ = "neil.millikin@gmail.com"
__copyright__ = "Copyright 2021, Neil Millikin"
__date__ = "2021/01/01"
__license__ = "GPLv3"
__version__ = "1.0.1"

import csv
import os
import sys
import inspect
import pprint

from itertools import islice
from collections import OrderedDict
from itertools import combinations

from PIL import Image, ImageDraw, ImageFont

from pixel_config import *

VERBOSITY = 2
font_library_path = "fonts"

helpful_debugging_utility_usages = """

       #  pretty print functions
       pp_2(str(inspect.stack()[0][2]), "Print something on onee line", here)

       # line print functions
       lp_2(str(inspect.stack()[0][2]), "print data structure, str(here))

   """
pp = pprint.PrettyPrinter(indent=4)

def lp_2(line_no, name, value):
    if VERBOSITY > 1:
        print("{0}_{1} = {2}\n".format(line_no, name, value))

def lp_3(line_no, name, value):
    if VERBOSITY > 2:
        print("{0}_{1} = {2}\n".format(line_no, name, value))


def pp_2(line_no, description, data_structure):
    if VERBOSITY > 1:
        print("\n{0} -- {1} ==>".format(line_no, description))
        pp.pprint(data_structure)
        print("")

def pp_3(line_no, description, data_structure):
    if VERBOSITY > 2:
        print("\n{0} -- {1} ==>".format(line_no, description))
        pp.pprint(data_structure)
        print("")

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))



CHROMOSOME_TO_RENDER = None

def get_match_pair_combinations(
        siblings_to_render,
        extra_match):

    all_matches_list = siblings_to_render

    match_pair_combinations = [
        (siblings_to_render[m[0]], siblings_to_render[m[1]]) for m in
        list(combinations(range(len(siblings_to_render)), 2))]
    if extra_match:
        for sib in siblings_to_render:
            match_pair_combinations.append((sib, extra_match))

        all_matches_list.append(extra_match)

    pp_2(str(inspect.stack()[0][2]), "match_pair_combinations",
         match_pair_combinations)

    return match_pair_combinations, all_matches_list


def get_match_pixel_dicts_for_siblings_to_render(
        all_matches_list):
    """
    For each match, read the raw file and construct a dictionary with this structure:

        { match_name_1 : { SNP_location_1 : ( tuple of SNP values A and B ),
                           SNP_location_2 : ( tuple of SNP values A and B ),
                           SNP_location_3 : ( tuple of SNP values A and B ),
                         },
          match_name_1 : { SNP_location_1 : ( tuple of SNP values A and B ),
                           SNP_location_2 : ( tuple of SNP values A and B ),
                           SNP_location_4 : ( tuple of SNP values A and B ),
                         },
         match_name_3 : { SNP_location_2 : ( tuple of SNP values A and B ),
                          SNP_location_3 : ( tuple of SNP values A and B ),
                          SNP_location_4 : ( tuple of SNP values A and B ),
                         },
        }

        Keep in mind that different testing companies include a different
        subset of SNP locations, so each match dict may have different keys.
    """
    match_SNP_values_dict = {}

    this_dir = os.path.dirname(os.path.realpath('__file__'))
    data_dir_name = DATA_FILE_DIRECTORY
    data_file_dir = os.path.join(this_dir, "{0}".format(data_dir_name))
    source_data_file_names = os.listdir(data_file_dir)

    raw_file_names = [f for f in source_data_file_names if 'raw' in f]

    for known_relative in all_matches_list:
        match_SNP_values_dict[known_relative] = {}
        try:
            kr_raw_file_name = [rfn for rfn in raw_file_names
                                if known_relative in rfn
                                and rfn.split('.')[ -1] == 'txt'][0]

            lp_2(str(inspect.stack()[0][2]), "kr_raw_file_name",
                 str(kr_raw_file_name))

        except IndexError as e:
            lp_2(str(inspect.stack()[0][2]),
                 "no such file found for {0}".format(known_relative), str(e))
            sys.exit()

        this_kr_raw_file = os.path.join(data_file_dir, kr_raw_file_name)

        with open(this_kr_raw_file, 'r', encoding='utf-8-sig') as raw_file:
            raw_reader = csv.reader([row for row in raw_file
                                     if not row.startswith ('#')
                                     and not row.startswith ('SNP')
                                     and not row.startswith('SNP')],
                                    delimiter='\t')
            for raw_row in raw_reader:

                if 'Ancestry' in kr_raw_file_name:
                    if str(raw_row[1]) != str(CHROMOSOME_TO_RENDER):
                        continue
                    raw_row_data = (raw_row[3], raw_row[4])

                else:
                    if str(raw_row[1]) != str(CHROMOSOME_TO_RENDER):
                        continue
                    raw_row_data = (raw_row[3][0], raw_row[3][1])

                match_SNP_values_dict[known_relative][
                    int(raw_row[2])] = raw_row_data

        pp_3(str(inspect.stack()[0][2]),
             "MATCH_PIXELS for {0}".format(known_relative),
             take(10, match_SNP_values_dict[known_relative].items()))

    return match_SNP_values_dict


def get_common_keys(
        match_SNP_values_dict):
    """
    Read and compare all matches, and get a list of ONLY those
    SNP locatiions ('rsid') that are represented in ALL of the matches.
    
    Otherwise you will get blank spots, chromosome images of different lenthts,
     or simply generate errors when trying to compare them.
    """
    list_of_match_pixel_dicts = []
    for kr_key in match_SNP_values_dict.keys():

        list_of_match_pixel_dicts.append(match_SNP_values_dict[kr_key])

    common_SNP_keys = set.intersection(*map(set, list_of_match_pixel_dicts))

    return common_SNP_keys


def get_common_key_SNP_dict(
        common_SNP_keys,
        match_SNP_values_dict):
    """
    Restructure the data dict to look like this, with only SNPs in all files:
    
    
        { SNP_location_1 : { match_name_1 : ( tuple of SNP values A and B ),
                             match_name_2 : ( tuple of SNP values A and B ),
                             match_name_3 : ( tuple of SNP values A and B ),
                           },
          SNP_location_2  : { match_name_1  : ( tuple of SNP values A and B ),
                             match_name_2 : ( tuple of SNP values A and B ),
                             match_name_3 : ( tuple of SNP values A and B ),
                           },
         SNP_location_3  : { match_name_1  : ( tuple of SNP values A and B ),
                             match_name_2 : ( tuple of SNP values A and B ),
                             match_name_3 : ( tuple of SNP values A and B ),
                           },
        }
    """
    common_key_SNP_dict = {}
    for common_SNP_key in common_SNP_keys:
        common_key_SNP_dict[common_SNP_key] = {}
        for kr_key in match_SNP_values_dict.keys():
            common_key_SNP_dict[common_SNP_key][kr_key] = \
                match_SNP_values_dict[kr_key][common_SNP_key]

    return common_key_SNP_dict


def insert_combo_match_type_into_common_key_SNP_dict(
        common_key_SNP_dict,
        match_pair_combinations):
    """
    transform the dataset into this structure, one entry for each SNP
    containing calculated match type (RED-YELLOW-GREEN):

    { SNP_1 : {
                'position': 289061,
                'JULIE': ('A', 'G'),
                'ALLISON': ('A', 'A'),
                'COLLETTE': ('A', 'G'),
                'J_A_Match': 'halfMatchSNP',
                'J_C_Match': 'fullMatchSNP',
                'A_C_Match': 'halfMatchSNP',
            },
     SNP_1 : {
                'position': 289061,
                'JULIE': ('G', 'G'),
                'ALLISON': ('T', 'T'),
                'COLLETTE': ('G', 'G'),
                'J_A_Match': 'noMatchSNP',
                'J_C_Match': 'fullMatchSNP',
                'A_C_Match': 'noMatchSNP',
            },
     SNP_1 : {
                'position': 289061,
                'JULIE': ('A', 'G'),
                'ALLISON': ('A', 'A'),
                'COLLETTE': ('A', 'G'),
                'J_A_Match': 'halfMatchSNP',
                'J_C_Match': 'fullMatchSNP',
                'A_C_Match': 'halfMatchSNP',
            },
    ),

    """

    useless_because_everything_is_identical_SNP_list = []
    """
    It is helpful to filter out 'noise' segments where all matches have the exact 
    same values.  These are SNPs where some variation is often found throughout 
    the human genome, but in this case there are no variants present in any match.

    There is no more reason to include these than for the millions of other SNP's 
    where variants either don't exist, are not commonly found, or are ubiquitious.
    
    The file has a setting for this, default is True
    """
    for SNP in common_key_SNP_dict:
        SNP_across_the_board_vals_list = []
        for sib in siblings_to_render:
            this_sib_vals = common_key_SNP_dict[SNP][sib][0:2]
            SNP_across_the_board_vals_list += this_sib_vals
        if len(set(SNP_across_the_board_vals_list)) == 1:
            useless_because_everything_is_identical_SNP_list.append(SNP)

    for SNP in common_key_SNP_dict:

        common_key_SNP_dict[SNP]['position'] = SNP

        for match_pair_combination in match_pair_combinations:
            mp_abbr = "{0}_{1}_Match".format(match_pair_combination[0][:3],
                                             match_pair_combination[1][:3])

            mpA = common_key_SNP_dict[SNP][
                "{0}".format(match_pair_combination[0])]
            mpB = common_key_SNP_dict[SNP][
                "{0}".format(match_pair_combination[1])]
            mpA_1 = mpA[0]
            mpA_2 = mpA[1]
            mpB_1 = mpB[0]
            mpB_2 = mpB[1]
            if (mpA_1 == mpB_1 and mpA_2 == mpB_2) or (
                    mpA_1 == mpB_2 and mpA_2 == mpB_1):
                match_type = 'fullMatchSNP'
            elif (mpA_1 == mpB_1 or mpA_2 == mpB_2) or (
                    mpA_1 == mpB_2 or mpA_2 == mpB_1):
                match_type = 'halfMatchSNP'
            else:
                match_type = 'noMatchSNP'

            common_key_SNP_dict[SNP][mp_abbr] = match_type

    if FILTER_COMPLETELY_MATCHED_SEGMENTS is True:

        no_across_the_board_matches_SNP_dict = {
                key: val for key, val
                in common_key_SNP_dict.items()
                if key not in useless_because_everything_is_identical_SNP_list}

        processed_and_sorted_SNP_dict = OrderedDict(
            sorted(no_across_the_board_matches_SNP_dict.items()))

    else:
        processed_and_sorted_SNP_dict = OrderedDict(
            sorted(common_key_SNP_dict.items()))

    pp_3(str(inspect.stack()[0][2]),
         "processed_and_sorted_SNP_dict",
         take(10, processed_and_sorted_SNP_dict.items()))

    return processed_and_sorted_SNP_dict


def create_comparison_base_strip_image(
        width,
        height):
    comarison_base_strip_image = Image.new(
            'RGB',
            (width, height),
            color='white')
    return comarison_base_strip_image


def create_chrom_whole_page_image(
        width,
        height,
        color):
    chromosome_full_page_image = Image.new(
            'RGB',
            (width, height),
            color=color)
    return chromosome_full_page_image


def draw_single_SNP_line(
        width,
        height,
        color):
    single_SNP_line = Image.new(
            'RGB',
            (width, height),
            color=color)
    return single_SNP_line


def show_match_graphics(
        processed_and_sorted_SNP_dict_table,
        match_pair_combinations):
    file_lines = len(processed_and_sorted_SNP_dict_table)

    if FILTER_COMPLETELY_MATCHED_SEGMENTS:
        title = "Pixel View Raw SNPs for Chr {0} -- filtered".format(
            CHROMOSOME_TO_RENDER)
    else:
        title = "Pixel View Raw SNPs for Chr {0} -- unfiltered".format(
            CHROMOSOME_TO_RENDER)

    matches_to_show = len(match_pair_combinations)
    chrom_whole_page_image = create_chrom_whole_page_image(

            max((file_lines
                 + CHROM_PAGE_LEFT_BORDER
                 + CHROM_PAGE_RIGHT_BORDER),
                MINIMUM_PAGE_WIDTH),

            (CHROM_PAGE_TOP_BORDER
            + CHROM_PAGE_BOTTOM_BORDER
            + CHROM_PAGE_TITLE_SPACE)
                + (matches_to_show * (
                    HEIGHT_OF_CHROMOSOME_IMAGE
                    + SPACE_BETWEEN_MATCHES)),

            BACKGROUND_COLOR)

    lp_3(str(inspect.stack()[0][2]), "matches_to_show", str(matches_to_show))

    page_draw = ImageDraw.Draw(chrom_whole_page_image)
    page_draw.text((CHROM_PAGE_LEFT_BORDER, CHROM_PAGE_TOP_BORDER), title,
            font=ImageFont.truetype(os.path.join(font_library_path, "Arial Bold.ttf"),
                                    CHROM_TITLE_TEXT_FONT_SIZE), fill='black')

    match_shown_number = 0

    for match_pair_combination in match_pair_combinations:

        mp_abbr = "{0}_{1}_Match".format(match_pair_combination[0][:3],
                                         match_pair_combination[1][:3])

        base_position = 0

        comparison_base_strip_image = create_comparison_base_strip_image(
                file_lines, HEIGHT_OF_CHROMOSOME_IMAGE)

        for SNP in processed_and_sorted_SNP_dict_table:

            relevant_match = processed_and_sorted_SNP_dict_table[SNP][mp_abbr]

            color = 'white'
            if relevant_match == 'noMatchSNP':
                color = NO_MATCH_SNP_COLOR
            elif relevant_match == 'halfMatchSNP':
                color = HALF_IDENTICAL_SNP_COLOR
            elif relevant_match == 'fullMatchSNP':
                color = FULLY_IDENTICAL_SNP_COLOR

            paste_position = (base_position, 0)

            single_pixel_line = draw_single_SNP_line(WIDTH_OF_SNP_LINE,
                    HEIGHT_OF_CHROMOSOME_IMAGE, color)

            comparison_base_strip_image.paste(single_pixel_line,
                    (paste_position))

            base_position += 1

        chrom_whole_page_image.paste(
                comparison_base_strip_image,

                (CHROM_PAGE_LEFT_BORDER,

                (CHROM_PAGE_TOP_BORDER + CHROM_PAGE_TITLE_SPACE
                    + (match_shown_number * (
                        HEIGHT_OF_CHROMOSOME_IMAGE
                        + SPACE_BETWEEN_MATCHES)))))

        draw = ImageDraw.Draw(chrom_whole_page_image)

        arial = ImageFont.truetype(
                os.path.join(font_library_path,
                             "Arial Bold.ttf"),
                CHROM_MATCH_TEXT_FONT_SIZE)

        draw.text((
                CHROM_PAGE_TEXT_BORDER,

                (CHROM_PAGE_TOP_BORDER
                    + CHROM_PAGE_TITLE_SPACE
                    + (match_shown_number
                       * (HEIGHT_OF_CHROMOSOME_IMAGE
                            + SPACE_BETWEEN_MATCHES)))),

                mp_abbr, font=arial, fill='black')

        match_shown_number += 1

    chrom_whole_page_image.show()


if __name__ == '__main__':

    def get_valid_chromosome_number():
        global CHROMOSOME_TO_RENDER
        SUGGESTED_CHROMOSOME_TO_RENDER = input("""
            "Enter a valid chromosome number to process 1-22  
            
            or type 'quit' to exit program     """)

        if SUGGESTED_CHROMOSOME_TO_RENDER == 'quit':
            print("goodbye")
            sys.exit()

        elif SUGGESTED_CHROMOSOME_TO_RENDER in [str(i) for i in range(1, 23)]:
            CHROMOSOME_TO_RENDER = int(SUGGESTED_CHROMOSOME_TO_RENDER)

        else:
            get_valid_chromosome_number()

    get_valid_chromosome_number()

    match_pair_combinations, all_matches_list \
        = get_match_pair_combinations(
            siblings_to_render,
            extra_match)

    match_SNP_values_dict = \
        get_match_pixel_dicts_for_siblings_to_render(
            all_matches_list)

    common_SNP_keys = \
        get_common_keys(
                match_SNP_values_dict)

    common_key_SNP_dict = \
        get_common_key_SNP_dict(
                common_SNP_keys,
                match_SNP_values_dict)

    processed_and_sorted_SNP_dict_table = \
        insert_combo_match_type_into_common_key_SNP_dict(
            common_key_SNP_dict,
            match_pair_combinations)

    show_match_graphics(
            processed_and_sorted_SNP_dict_table,
            match_pair_combinations)
