//
// Created by Zhenhao on 8/5/2022.
//

/**
 * This file is adopted from SeqAn3 tutorial: https://docs.seqan.de/seqan/3.0.2/tutorial_read_mapper.html.
 */

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <span>

#include "../hurdle_matrix.h"
#include "../seqan3_main.h"

struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<seqan3::dna5>> seqs;
};

void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & storage)
{
    seqan3::sequence_file_input reference_in{reference_path};
    for (auto & [seq, id, qual] : reference_in)
    {
        storage.ids.push_back(std::move(id));
        storage.seqs.push_back(std::move(seq));
    }
}

void map_reads(std::filesystem::path const & query_path,
               std::filesystem::path const & index_path,
               std::filesystem::path const & sam_path,
               reference_storage_t & storage,
               uint8_t const errors)
{
    // we need the alphabet and text layout before loading
    seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> index;
    {
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }

    seqan3::sequence_file_input query_file_in{query_path};

    seqan3::sam_file_output sam_out{sam_path, seqan3::fields<seqan3::field::seq,
            seqan3::field::id,
            seqan3::field::ref_id,
            seqan3::field::ref_offset,
            seqan3::field::alignment,
            seqan3::field::qual,
            seqan3::field::mapq>{}};

    seqan3::configuration const search_config = seqan3::search_cfg::max_error_total{
            seqan3::search_cfg::error_count{errors}} |
                                                seqan3::search_cfg::hit_single_best{};

    hurdle_matrix<int_128bit>* matrix = new hurdle_matrix<int_128bit>(GLOBAL, 1, 1, 1);

    for (auto && record : query_file_in)
    {
        auto & query = record.sequence();
        //std::cout << std::string(query) << std::endl;
        for (auto && result : search(query, index, search_config))
        {
            size_t start = result.reference_begin_position() ? result.reference_begin_position() - 1 : 0;
            std::span text_view{std::data(storage.seqs[result.reference_id()]) + start, query.size() + 1};
            // Run the hurdle matrix
            matrix->reset(seqan_dna_to_cstring(query).c_str(),
                          query.size(),
                          seqan_dna_to_cstring(text_view).c_str(),
                          text_view.size(), 3);
            matrix->run();
            std::pair<std::span<seqan3::dna5>, std::span<seqan3::dna5>> alignment = {query, text_view};

            sam_out.emplace_back(query,
                                 record.id(),
                                 storage.ids[result.reference_id()],
                                 start + 2, // TODO: correct this number
                                 alignment, // FIXME: let hurdle matrix return a comparable type to seqan3::gap
                                 record.base_qualities(),
                                 60u + matrix->get_cost()
                                 );
        }
    }
}

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & index_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    reference_storage_t storage{};
    read_reference(reference_path, storage);
    map_reads(query_path, index_path, sam_path, storage, errors);
}

struct cmd_arguments
{
    std::filesystem::path reference_path{};
    std::filesystem::path query_path{};
    std::filesystem::path index_path{};
    std::filesystem::path sam_path{"out.sam"};
    uint8_t errors{0};
};

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli";
    parser.info.short_description = "Map reads against a reference.";
    parser.info.version = "1.0.0";
    parser.add_option(args.reference_path, 'r', "reference", "The path to the reference.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa","fasta"}});
    parser.add_option(args.query_path, 'q', "query", "The path to the query.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fq","fastq"}});
    parser.add_option(args.index_path, 'i', "index", "The path to the index.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"index", "fai"}});
    parser.add_option(args.sam_path, 'o', "output", "The output SAM file path.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sam"}});
    parser.add_option(args.errors, 'e', "error", "Maximum allowed errors.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 4});
}

int main(int argc, char const ** argv)
{
    seqan3::argument_parser parser("Mapper", argc, argv);
    cmd_arguments args{};

    initialise_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    run_program(args.reference_path, args.query_path, args.index_path, args.sam_path, args.errors);

    return 0;
}