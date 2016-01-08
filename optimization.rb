#!/usr/bin/env ruby

# This library is used to generate optimization problem
#
# Author:: Matteo Ragni
# Copyright:: Copyright (2015) University of Trento
# License:: Distributes under the same terms as Ruby

gems = %w|nmatrix nmatrix/lapacke colorize pry-byebug|
gems.each do |gem|
  #raise LoadError, "Cannot Load Gem :: #{gem}" unless require gem
  require gem
end

##
# Full container for a set of optimization problem completely defined i Ruby.
# Even if it is not convenient in term of computational time. This is for
# sure a very confortable way to prototype algorigthms.
module Optimization
  load './functions.rb'

  ##
  # This class is a general container for the algorithm that I will define
  class Algorithm
    # Objective function
    attr_reader :objective
    # List of constraints for the problem
    attr_reader :constraints
    # History of the optimization problem
    attr_reader :history
    # Options for the algorithm
    attr_accessor :options
    # Size for problem
    attr_reader :size
    # Number of equalities and inequalities
    attr_reader :eq_size, :ineq_size

    ##
    # Initialize a new optimization problem. As for now, multiple dimensions
    # linear functions are not supported
    def initialize(options, objective)
      raise ArgumentError, "options must be an Hash" unless options.is_a? Hash
      raise ArgumentError, "Objective function must be of type ObjectiveFunction" unless
        (objective.is_a? ObjectiveFunction or
         objective.is_a? LinearObjectiveFunction or
         objective.is_a? QuadraticObjectiveFunction)

      @options = options
      @options[:debug]     = 0     unless @options[:debug]
      @options[:tolerance] = 1E-10 unless @options[:tolerance]

      @objective = objective
      @size = objective.size
      @constraints = []
      @eq_size = 0
      @ineq_size = 0
    end

    ##
    # Add a new constraint to the problem
    def add_constraint(c)
      raise ArgumentError, "Argument must be a constraint function" unless
        (c.is_a? ConstraintFunction or
         c.is_a? LinearConstraintFunction or
         c.is_a? QuadraticConstraintFunction)
      raise ArgumentError, "Function must accept x of size #{@size}. Is of size #{c.size}" unless @size == c.size
      @eq_size += 1 if c.type == :equality
      @ineq_size += 1 if c.type == :inequality
      @constraints << c
    end

    ##
    # Clear all the constraints resetting counters
    def remove_constraints
      @constraints = []
      @eq_size     = 0
      @ineq_size   = 0
    end

    ##
    # Solves the optimization problem given a starting guess
    def solve(x0)
      raise RuntimeError, "This method is not implemented!"
    end

    ##
    # Returns false if a point is not feasible (does not live in _P_ domain)
    def feasible?(x)
      feas = true
      @constraints.each do |c|
        feas = false if is_violated?(c,x)
      end
      #puts inspect_constraints(x) if @options[:debug] > 4
      info(inspect_constraints(x), 4, __method__, __LINE__)
      return feas
    end

    ##
    # Returns true if a point is infeasible
    def infeasible?(x)
      cnt = []
      @constraints.each do |c|
        cnt << is_violated?(c,x)
      end
      puts inspect_constraints(x) if @options[:debug] > 4
      return cnt.include? true
    end

    ##
    # Evaluated if a constraints is violated with respect to
    # configured error
    def is_violated?(c, x)
      if c.type == :equality
        return ((c.f(x)).abs > @options[:tolerance])
      else
        return (c.f(x) < -@options[:tolerance])
      end
    end

    ##
    # Inspect constraint status
    def inspect_constraints(x = nil)
      output = "Constraints status\n"
      @constraints.each_with_index do |c, i|
        output += " * #{i+1}) c(#{x ? x.to_flat_array : 'x'}) #{c.type == :equality ? '=' : '≥'} 0? #{x ? (is_violated?(c,x) ? 'violated' : 'not violated') : ''} " + "[ c(x) = #{c.f(x)}]" + "\n"
      end
      return output
    end

    def info(msg, level = 1, method = nil, line = nil)
      color = case level
      when 0
       :light_yellow
      when :result
       :light_yellow
       level = 0
      when 1
       :light_white
      when 2
       :light_blue
      when 3
       :light_green
      when 4
       :yellow
      when 5
       :light_red
      else
        :white
      end
      puts "#{method ? "#{method}:" : ""}#{line ? "#{line}:" : ""} " +
            (level == 0 ? " --> ".colorize(:light_white) : "") +
            msg.to_s.colorize(color) if @options[:debug] > level
    end
  end

  load './general_quadratic_optimizer.rb'
end

if $0 == __FILE__ then

  # Objective function
  m = NMatrix.new [2,2], [3.0, -1.0, -2.0, 7.0]
  g = NMatrix.new [2,1], [2.0, -3.0]
  c = 4.0
  obj = Optimization::QuadraticObjectiveFunction.new m, g, c

  # Constraint :: y + 1 >= 0
  a1 = NMatrix.new [2,1], [0.0, 1.0]
  b1 = 1.0
  cnt1 = Optimization::LinearConstraintFunction.new :inequality, a1, b1

  # Constraint :: x + 1/4 >= 0
  a2 = NMatrix.new [2,1], [1.0, 0.0]
  b2 = 0.25
  cnt2 = Optimization::LinearConstraintFunction.new :inequality, a2, b2

  # Constraint :: - x - y >= 0
  a3 = NMatrix.new [2,1], [0.0, 1.0]
  b3 = 0.0
  cnt3 = Optimization::LinearConstraintFunction.new :inequality, a3, b3

  # Constraint :: 12/19 * x + y = 0
  a4 = NMatrix.new [2,1], [12.0/19.0, 1.0]
  b4 = 0.0
  cnt4 = Optimization::LinearConstraintFunction.new :equality, a4, b4

  # Creation of a new test scenario
  alg = Optimization::GeneralQuadraticOptimizer.new({debug: 10, iterations: 100}, obj)
  alg.add_constraint(cnt1)
  alg.add_constraint(cnt2)
  alg.add_constraint(cnt3)
  #alg.add_constraint(cnt4)

  x0 = NMatrix.new [2,1], [0.0, 0.0]
  x = alg.solve x0

=begin
  Q = NMatrix.new [3,3],
    [ 1.0, 0.0, 0.0,
      0.0, 2.0, 0.0,
      0.0, 0.0, 3.0]
  G = NMatrix.new [3,1], [ 3.0, 2.0, 1.0 ]
  C = 1.0
  obj = Optimization::QuadraticObjectiveFunction.new Q, G, C

  A1 = NMatrix.new [3,1], [ 1.0, 2.0, 3.0 ]
  B1 = 0.0
  c1 = Optimization::LinearConstraintFunction.new :equality, A1, B1

  alg = Optimization::QuadraticOptimizer.new({debug: 10}, obj)
  alg.add_constraint c1

  x, lambdas = alg.solve
  puts " -> x = #{x.to_flat_array}".yellow
=end
  #binding.pry
end
