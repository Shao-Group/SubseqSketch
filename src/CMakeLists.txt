add_library(fasta_reader fasta_reader.cpp)

add_library(subsequences subsequences.cpp)
target_link_libraries(subsequences PUBLIC fasta_reader)

# add_library(tokenized_sequence tokenized_sequence.cpp)

add_library(sss_array sss_array.cpp)
target_link_libraries(sss_array PUBLIC OpenMP::OpenMP_CXX)

add_executable(SubseqSketch subseq_sketch.cpp)
target_link_libraries(SubseqSketch PRIVATE subsequences)
# target_link_libraries(SubseqSketch PRIVATE tokenized_sequence)
target_link_libraries(SubseqSketch PRIVATE sss_array)

