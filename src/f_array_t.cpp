#include "f_array_t.hpp"

f_array_t::size_type f_array_t::serialize(std::ostream &out, structure_tree_node *v, std::string name) const{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_C.serialize(out, child, "C");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void f_array_t::load(std::istream &in){
    m_C.load(in);
}

void f_array_t::swap(f_array_t &in){
    m_C.swap(in.m_C);
}


