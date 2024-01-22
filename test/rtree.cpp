#include <gtest/gtest.h>

#include <RTree.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <memory>

namespace er = eh::rtree;

TEST( RTreeTest, Insert )
{
  using rtree_type = er::RTree<er::aabb_t<int>,er::aabb_t<int>,int>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;


  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -1000, 1000 );
  rtree_type rtree;

  for( int i=0; i<1000; ++i )
  {
    int min_ = dist(mt);
    int max_ = dist(mt);
    if( max_ < min_ ){ std::swap(min_,max_); }
    rtree.insert( {{min_,max_}, i} );
    ASSERT_EQ( rtree.size(), i+1 );

    // entries count check
    for( int level=0; level<rtree.leaf_level(); ++level )
    {
      for( auto ni=rtree.begin(level); ni!=rtree.end(level); ++ni )
      {
        if( level != 0 )
        {
          ASSERT_GE( ni->size(), rtree_type::MIN_ENTRIES ) << "level " << level << " below MIN_ENTRIES";
        }
        ASSERT_LE( ni->size(), rtree_type::MAX_ENTRIES ) << "level " << level << " above MAX_ENTRIES";
      }
    }
    for( auto ni=rtree.leaf_begin(); ni!=rtree.leaf_end(); ++ni )
    {
      if( rtree.leaf_level() != 0 )
      {
        ASSERT_GE( ni->size(), rtree_type::MIN_ENTRIES ) << "leaf below MIN_ENTRIES";
      }
      ASSERT_LE( ni->size(), rtree_type::MAX_ENTRIES ) << "leaf above MAX_ENTRIES";
    }
    // entries bound check
    for( int level=0; level<rtree.leaf_level(); ++level )
    {
      for( auto ni=rtree.begin(level); ni!=rtree.end(level); ++ni )
      {
        auto *node = *ni;
        for( auto &c : *node )
        {
          ASSERT_EQ( c.second->parent(), node );
          auto a = c.first;
          if( level+1==rtree.leaf_level() )
          {
            auto b = c.second->as_leaf()->calculate_bound();
            ASSERT_TRUE( traits::is_inside(a,b) ) << "i: " << i << ", level: " << level;
          }else {
            auto b = c.second->as_node()->calculate_bound();
            ASSERT_TRUE( traits::is_inside(a,b) ) << "i: " << i << ", level: " << level;
          }
        }
      }
    }

    // all data valid check
    std::vector<bool> valid( i+1, false );
    for( auto x : rtree )
    {
      ASSERT_FALSE( valid[x.second] ) << "i: " << i << ", " << "x: " << x.second << " already exist";
      valid[x.second] = true;
    }
    for( int x=0; x<=i; ++x )
    {
      ASSERT_TRUE( valid[x] ) << "i: " << i << ", " << "x: " << x << " not exist";
    }
  }
}

TEST( RTreeTest, Erase )
{
  using rtree_type = er::RTree<er::aabb_t<int>,er::aabb_t<int>,int>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;


  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -1000, 1000 );
  rtree_type rtree;

  for( int i=0; i<1000; ++i )
  {
    int min_ = dist(mt);
    int max_ = dist(mt);
    if( max_ < min_ ){ std::swap(min_,max_); }
    rtree.insert( {{min_,max_}, i} );
  }

  std::vector<bool> data_inserted( 1000, true );
  for( int i=0; i<1000; ++i )
  {
    const int cur_size = 1000-i;
    int erase_index = std::uniform_int_distribution<int>( 0, cur_size-1 )( mt );
    auto it = rtree.begin();
    std::advance( it, erase_index );
    ASSERT_GE( it->second, 0 );
    ASSERT_LE( it->second, 1000 );
    ASSERT_TRUE( data_inserted[it->second] ) << it->second << ": already deleted";
    data_inserted[it->second] = false;
    rtree.erase( it );
    ASSERT_EQ( rtree.size(), cur_size-1 );

    // entries count check
    for( int level=0; level<rtree.leaf_level(); ++level )
    {
      for( auto ni=rtree.begin(level); ni!=rtree.end(level); ++ni )
      {
        if( level != 0 )
        {
          ASSERT_GE( ni->size(), rtree_type::MIN_ENTRIES ) << "level " << level << " below MIN_ENTRIES";
        }
        ASSERT_LE( ni->size(), rtree_type::MAX_ENTRIES ) << "level " << level << " above MAX_ENTRIES";
      }
    }
    for( auto ni=rtree.leaf_begin(); ni!=rtree.leaf_end(); ++ni )
    {
      if( rtree.leaf_level() != 0 )
      {
        ASSERT_GE( ni->size(), rtree_type::MIN_ENTRIES ) << "leaf below MIN_ENTRIES";
      }
      ASSERT_LE( ni->size(), rtree_type::MAX_ENTRIES ) << "leaf above MAX_ENTRIES";
    }

    // entries bound check
    for( int level=0; level<rtree.leaf_level(); ++level )
    {
      for( auto ni=rtree.begin(level); ni!=rtree.end(level); ++ni )
      {
        auto *node = *ni;
        for( auto &c : *node )
        {
          ASSERT_EQ( c.second->parent(), node );
          auto a = c.first;
          if( level+1==rtree.leaf_level() )
          {
            auto b = c.second->as_leaf()->calculate_bound();
            ASSERT_TRUE( traits::is_inside(a,b) ) << "i: " << i << ", level: " << level;
          }else {
            auto b = c.second->as_node()->calculate_bound();
            ASSERT_TRUE( traits::is_inside(a,b) ) << "i: " << i << ", level: " << level;
          }
        }
      }
    }
    ASSERT_EQ( std::distance(rtree.begin(),rtree.end()), cur_size-1 ) << "distance(tree.iterator) is not equal to " << i+1;
  }
  for( int i=0; i<1000; ++i )
  {
    ASSERT_FALSE( data_inserted[i] ) << "data " << i << " exist after deletion";
  }
}

TEST( RTreeTest, Search )
{
  using rtree_type = er::RTree<er::aabb_t<int>,er::aabb_t<int>,int>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;
  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -100, 100 );

  rtree_type rtree;
  bound_type search_range = { -10, 10 };

  std::vector<int> inside_list, overlap_list;

  for( int i=0; i<1000; ++i )
  {
    int min_ = dist(mt);
    int max_ = dist(mt);
    if( max_ < min_ ){ std::swap(min_,max_); }
    rtree.insert( {{min_,max_}, i} );

    if( traits::is_inside(search_range, bound_type{min_,max_}) )
    {
      inside_list.push_back( i );
    }
    if( traits::is_overlap(search_range, bound_type{min_,max_}) )
    {
      overlap_list.push_back( i );
    }

    std::vector<int> cur_inside, cur_overlap;
    rtree.search_inside( search_range,
      [&]( rtree_type::value_type v )
      {
        cur_inside.push_back( v.second );
        return false;
      }
    );
    std::sort( cur_inside.begin(), cur_inside.end() );
    ASSERT_EQ( inside_list.size(), cur_inside.size() );
    for( int j=0; j<cur_inside.size(); ++j )
    {
      ASSERT_EQ( cur_inside[j], inside_list[j] );
    }

    rtree.search_overlap( search_range,
      [&]( rtree_type::value_type v )
      {
        cur_overlap.push_back( v.second );
        return false;
      }
    );

    std::sort( cur_overlap.begin(), cur_overlap.end() );
    ASSERT_EQ( overlap_list.size(), cur_overlap.size() );
    for( int j=0; j<cur_overlap.size(); ++j )
    {
      ASSERT_EQ( cur_overlap[j], overlap_list[j] );
    }
  }
}

TEST( RTreeTest, Assign )
{
  using rtree_type = er::RTree<er::aabb_t<int>,int,int>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;
  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -100, 100 );

  rtree_type rtree;
  std::vector<rtree_type::value_type> original;

  for( int i=0; i<1000; ++i )
  {
    int pos = dist(mt);
    rtree.insert( {pos, i} );
    original.push_back( {pos,i} );
  }
  std::sort( original.begin(), original.end() );

  rtree_type rtree_copy = rtree;
  std::vector<rtree_type::value_type> copied( rtree_copy.begin(), rtree_copy.end() );
  std::sort( copied.begin(), copied.end() );

  ASSERT_EQ( original.size(), copied.size() );
  for( int i=0; i<original.size(); ++i )
  {
    ASSERT_EQ( original[i], copied[i] );
  }

  rtree_type rtree_moved = std::move(rtree);
  std::vector<rtree_type::value_type> moved( rtree_moved.begin(), rtree_moved.end() );
  std::sort( moved.begin(), moved.end() );
  ASSERT_EQ( original.size(), copied.size() );
  for( int i=0; i<original.size(); ++i )
  {
    ASSERT_EQ( original[i], moved[i] );
  }
}

TEST( RTreeTest, Clear )
{
  using rtree_type = er::RTree<er::aabb_t<int>,int,int>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;
  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -100, 100 );

  rtree_type rtree;
  std::vector<rtree_type::value_type> original;

  for( int i=0; i<1000; ++i )
  {
    int pos = dist(mt);
    rtree.insert( {pos, i} );
    original.push_back( {pos,i} );
  }
  std::sort( original.begin(), original.end() );

  rtree.clear();
  ASSERT_EQ( rtree.size(), 0 );
  ASSERT_EQ( rtree.leaf_level(), 0 );
  ASSERT_EQ( rtree.begin(), rtree.end() );
}

TEST( RTreeTest, MoveOnlyValue )
{
  std::vector<int> destroyed( 1000, 0 );
  struct move_only_t
  {
    int i;
    std::vector<int> *destroyed = nullptr;

    void destroy()
    {
      if( i != -1 )
      {
        (*destroyed)[i]++;
      }
    }
    move_only_t( int i_, std::vector<int>&d ) : i(i_), destroyed(&d) {}
    move_only_t( move_only_t && rhs )
    {
      i = rhs.i;
      destroyed = rhs.destroyed;
      rhs.i = -1;
    }
    move_only_t( const move_only_t & ) = delete;

    ~move_only_t()
    {
      destroy();
    }

    move_only_t& operator=( move_only_t && rhs )
    {
      destroy();
      i = rhs.i;
      destroyed = rhs.destroyed;
      rhs.i = -1;
      return *this;
    }
    move_only_t& operator=( const move_only_t & ) = delete;
  };
  using rtree_type = er::RTree<er::aabb_t<int>,int,move_only_t>;
  using traits = rtree_type::traits;
  using bound_type = rtree_type::geometry_type;
  using node_type = rtree_type::node_type;


  std::mt19937 mt( std::random_device{}() );
  std::uniform_int_distribution<int> dist( -1000, 1000 );
  {
    rtree_type rtree;

    for( int i=0; i<1000; ++i )
    {
      int p = dist(mt);
      rtree.insert( {p, move_only_t(i,destroyed)} );
      ASSERT_EQ( rtree.size(), i+1 );

      // all data valid check
      std::vector<bool> valid( i+1, false );
      for( auto &x : rtree )
      {
        ASSERT_FALSE( valid[x.second.i] ) << "i: " << i << ", " << "x: " << x.second.i << " already exist";
        valid[x.second.i] = true;
      }
      for( int x=0; x<=i; ++x )
      {
        ASSERT_TRUE( valid[x] ) << "i: " << i << ", " << "x: " << x << " not exist";
      }
    }
    // destruction check
  }
  for( int i=0; i<1000; ++i )
  {
    ASSERT_EQ( destroyed[i], 1 ) << "i: " << i << " destroyed " << destroyed[i] << " times";
  }
}