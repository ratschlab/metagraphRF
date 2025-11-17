#include "metagraph.hpp"
#include "cli/load/load_annotation.hpp"
#include "seq_io/sequence_io.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "prelude.h"
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/io.h"

namespace py = pybind11;
using namespace std;

struct Alignment {
    string ctg;
    int r_st = 0, r_en = 0, strand = 1;
    float pres_frac = 0.0f;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

struct Request {
    int channel = 0;
    string id;
    string seq;
    Request() = default;
    Request(int channel, string &id, string &seq): channel(channel), id(id), seq(seq) {}
};

struct Response {
    int channel = 0;
    string id;
    Alignment alignment;
    Response() = default;
    Response(int channel, string &id, Alignment &alignment): channel(channel), id(id), alignment(alignment) {}
};

struct ResponseGenerator {
    parlay::sequence<Response> responses;
    explicit ResponseGenerator(parlay::sequence<Response> &responses): responses(responses) {}
    Response next() {
        if (responses.empty())
            throw py::stop_iteration();
        else {
            auto response = responses.back();
            responses.pop_back();
            return response;
        }
    }
    ResponseGenerator& iter() {
        return *this;
    }
};

struct Index {
    Index(const string& graph_path, const string& annotation_path, double pres_frac = 0.3):
    pres_frac(pres_frac) {
        log_info("Loading graph..");
        graph = mtg::cli::load_critical_dbg(graph_path);
        if (!graph) log_error("Failed to load graph from %s", graph_path.c_str());
        log_info("Loading annotations...");

        auto anno_type = mtg::cli::parse_annotation_type(annotation_path);
        annotation = mtg::cli::initialize_annotation(anno_type);
        if (!annotation->load(annotation_path)) log_error("Failed to load annotation from %s", annotation_path.c_str());

        anno_graph = make_unique<mtg::graph::AnnotatedDBG>(graph, move(annotation));
    }

    Alignment query(string &sequence) {
        vector<mtg::graph::DeBruijnGraph::node_index> mapping;
        anno_graph->get_graph().map_to_nodes_sequentially(sequence, [&](auto node) {
            mapping.push_back(node);
        });
        // Count how many k-mers matched
        size_t total_kmers = mapping.size();
        size_t matched_kmers = 0;
        for (const auto& node : mapping) {
            if (node != mtg::graph::DeBruijnGraph::npos) {
                matched_kmers++;
            }
        }

        double match_fraction = total_kmers > 0
            ? static_cast<double>(matched_kmers) / total_kmers
            : 0.0;

        if (match_fraction >= pres_frac && matched_kmers > 0) {
            // Get annotations for matched nodes
            map<string, size_t> label_counts;

            for (const auto& node : mapping) {
                if (node != mtg::graph::DeBruijnGraph::npos) {
                    // Convert graph node index to annotation index
                    auto anno_index = mtg::graph::AnnotatedDBG::graph_to_anno_index(node);
                    auto labels = anno_graph->get_annotator().get_labels(anno_index);
                    for (const auto& label : labels) {
                        label_counts[label]++;
                    }
                }
            }

            if (!label_counts.empty()) {
                using pair_type = decltype(label_counts)::value_type;
                auto best_label_it = max_element(label_counts.begin(), label_counts.end(),
                    [](const pair_type& p1, const pair_type& p2) { return p1.second < p2.second; });
                return {best_label_it->first.c_str(), true, 0,
                    static_cast<float>(match_fraction), static_cast<int>(sequence.length())};
            }
            return {"*", true, 0, 0, static_cast<int>(sequence.length())};
        }
        return {"*", true, 0, 0, static_cast<int>(sequence.length())};
    }

    vector<Alignment> query_batch(const py::list& sequences) {
        auto nr = sequences.size();
        vector<Alignment> results(nr);
        parlay::for_each(parlay::iota(nr), [&](size_t i){
            auto sequence = sequences[i].cast<string>();
            results[i] = query(sequence);
        });
        return results;
    }

    ResponseGenerator query_stream(const py::iterator& reads) {
        parlay::sequence<Request> requests;
        for (auto &read: reads) {
            auto request = read.cast<Request>();
            requests.push_back(request);
        }
        auto responses = parlay::tabulate(requests.size(), [&](size_t i) {
            auto alignment = query(requests[i].seq);
            return Response(requests[i].channel, requests[i].id, alignment);
        });
        return ResponseGenerator(responses);
    }

    shared_ptr<mtg::graph::DeBruijnGraph> graph;
    unique_ptr<mtg::annot::MultiLabelAnnotation<string>> annotation;
    unique_ptr<mtg::graph::AnnotatedDBG> anno_graph;
    double pres_frac;

};

PYBIND11_MODULE(_core, m) {
    py::class_<Alignment>(m, "Alignment")
            .def(py::init<>())  // Default constructor
            .def(py::init<const char*, bool, int, float, int>(),  // Parameterized constructor
                 py::arg("header"), py::arg("fwd"), py::arg("start"),
                 py::arg("pres_frac"), py::arg("qry_len"))
            .def_readonly("ctg", &Alignment::ctg)
            .def_readonly("r_st", &Alignment::r_st)
            .def_readonly("r_en", &Alignment::r_en)
            .def_readonly("strand", &Alignment::strand)
            .def_readonly("pres_frac", &Alignment::pres_frac);

    py::class_<Index>(m, "Index")
            .def(py::init<const string&, const string&, double>(),
                py::arg("graph"), py::arg("annotation"),
                py::arg("pres_frac") = 0.3)
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch)
            .def("query_stream", &Index::query_stream);

    py::class_<Request>(m, "Request")
            .def(py::init<int, string&, string&>(), py::arg("channel"), py::arg("id"), py::arg("seq"))
            .def_readwrite("channel", &Request::channel)
            .def_readwrite("id", &Request::id)
            .def_readwrite("seq", &Request::seq);

    py::class_<Response>(m, "Response")
            .def(py::init<int, string&, Alignment&>(), py::arg("channel"), py::arg("id"), py::arg("alignment"))
            .def_readwrite("channel", &Response::channel)
            .def_readwrite("id", &Response::id)
            .def_readwrite("alignment", &Response::alignment);

    py::class_<ResponseGenerator>(m, "ResponseGenerator")
            .def("__iter__", &ResponseGenerator::iter)
            .def("__next__", &ResponseGenerator::next);
}