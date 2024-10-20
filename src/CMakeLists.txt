include_directories(${SubseqEmbed_SOURCE_DIR}/include)
add_library(fasta_reader fasta_reader.cpp)

add_library(subsequences subsequences.cpp)
target_link_libraries(subsequences PUBLIC fasta_reader)

add_library(tokenized_sequence tokenized_sequence.cpp)

add_library(rssebd_array rssebd_array.cpp)

add_executable(SubseqEmbed subseq_embed.cpp)
target_link_libraries(SubseqEmbed PRIVATE subsequences)
target_link_libraries(SubseqEmbed PRIVATE tokenized_sequence)
target_link_libraries(SubseqEmbed PRIVATE rssebd_array)

