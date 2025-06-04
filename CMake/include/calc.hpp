#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <expected>
#include <stack>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace calc {

enum class ArithmeticError {
  Overflow,
  Underflow,
  InvalidInput,
  OutOfRange,
  DivisionByZero,
  InvalidExpression,
  UnbalancedParentheses,
  UnknownOperator,
  DomainError
};

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <Arithmetic T>
std::expected<T, ArithmeticError> parse_number(std::string_view input) {
  auto trimmed = input;
  trimmed.remove_prefix(
      std::min(trimmed.find_first_not_of(" \t"), trimmed.size()));
  trimmed.remove_suffix(trimmed.size() - trimmed.find_last_not_of(" \t") - 1);

  if (trimmed.empty()) {
    return std::unexpected(ArithmeticError::InvalidInput);
  }

  std::string normalized(trimmed);
  std::replace(normalized.begin(), normalized.end(), ',', '.');

  T value;
  auto [ptr, ec] = std::from_chars(
      normalized.data(), normalized.data() + normalized.size(), value);

  if (ec == std::errc::invalid_argument) {
    return std::unexpected(ArithmeticError::InvalidInput);
  }
  if (ec == std::errc::result_out_of_range) {
    return std::unexpected(ArithmeticError::OutOfRange);
  }
  if (ptr != normalized.data() + normalized.size()) {
    return std::unexpected(ArithmeticError::InvalidInput);
  }

  return value;
}

namespace detail {

constexpr double PI = 3.14159265358979323846;
constexpr double E = 2.71828182845904523536;

constexpr double abs(double x) { return x >= 0 ? x : -x; }

constexpr double factorial(int n) {
  if (n < 0)
    return 0;
  double result = 1;
  for (int i = 2; i <= n; ++i)
    result *= i;
  return result;
}

inline double sin_approx(double x) {
  x = std::fmod(x, 2 * PI);
  if (x < 0)
    x += 2 * PI;

  double result = 0.0;
  double term = x;
  double x2 = x * x;
  int n = 1;
  while (abs(term) > 1e-10 && n < 20) {
    result += term;
    n += 2;
    term *= -x2 / (n * (n - 1));
  }
  return result;
}

constexpr double cos_approx(double x) {
  while (x > PI)
    x -= 2 * PI;
  while (x < -PI)
    x += 2 * PI;

  double result = 0.0;
  double term = 1.0;
  double x2 = x * x;
  int n = 0;
  while (abs(term) > 1e-10 && n < 20) {
    result += term;
    n += 2;
    term *= -x2 / (n * (n - 1));
  }
  return result;
}

constexpr std::expected<double, ArithmeticError> tan_approx(double x) {
  double cos_val = cos_approx(x);
  if (abs(cos_val) < 1e-10) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return sin_approx(x) / cos_val;
}

constexpr std::expected<double, ArithmeticError> cot_approx(double x) {
  auto tan_val = tan_approx(x);
  if (!tan_val || abs(*tan_val) < 1e-10) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return 1.0 / *tan_val;
}

constexpr std::expected<double, ArithmeticError> sec_approx(double x) {
  double cos_val = cos_approx(x);
  if (abs(cos_val) < 1e-10) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return 1.0 / cos_val;
}

inline std::expected<double, ArithmeticError> csc_approx(double x) {
  const double sin_val = sin_approx(x);
  if (abs(sin_val) < 1e-10) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return 1.0 / sin_val;
}

constexpr std::expected<double, ArithmeticError> asin_approx(double x) {
  if (x < -1.0 || x > 1.0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  double result = 0.0;
  double term = x;
  double x2 = x * x;
  int n = 1;
  while (abs(term) > 1e-10 && n < 20) {
    result += term;
    n += 2;
    term *= (x2 * (n - 2)) / n * (n - 1) / (n + 1);
  }
  return result;
}

constexpr std::expected<double, ArithmeticError> acos_approx(double x) {
  auto asin_val = asin_approx(x);
  if (!asin_val)
    return asin_val;
  return PI / 2.0 - *asin_val;
}

constexpr double atan_approx(double x) {
  if (abs(x) > 1.0) {
    double inner = atan_approx(1.0 / x);
    return x > 0 ? PI / 2.0 - inner : -PI / 2.0 - inner;
  }
  double result = 0.0;
  double term = x;
  double x2 = x * x;
  int n = 1;
  while (abs(term) > 1e-10 && n < 20) {
    result += term;
    n += 2;
    term *= -x2 * (n - 2) / n;
  }
  return result;
}

constexpr std::expected<double, ArithmeticError> acot_approx(double x) {
  if (x == 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return PI / 2.0 - atan_approx(x);
}

constexpr std::expected<double, ArithmeticError> asec_approx(double x) {
  if (x > -1.0 && x < 1.0 && x != 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return acos_approx(1.0 / x);
}

constexpr std::expected<double, ArithmeticError> acsc_approx(double x) {
  if (x > -1.0 && x < 1.0 && x != 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  return asin_approx(1.0 / x);
}

constexpr std::expected<double, ArithmeticError> log10_approx(double x) {
  if (x <= 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }

  double z = (x - 1.0) / (x + 1.0);
  double z2 = z * z;
  double result = 0.0;
  double term = z;
  int n = 1;
  while (abs(term) > 1e-10 && n < 50) {
    result += term / n;
    n += 2;
    term *= z2 * (n - 2);
  }
  result *= 2.0;
  constexpr double ln10 = 2.30258509299404568402;
  return result / ln10;
}

constexpr std::expected<double, ArithmeticError> ln_approx(double x) {
  if (x <= 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  double z = (x - 1.0) / (x + 1.0);
  double z2 = z * z;
  double result = 0.0;
  double term = z;
  int n = 1;
  while (abs(term) > 1e-10 && n < 50) {
    result += term / n;
    n += 2;
    term *= z2 * (n - 2);
  }
  return 2.0 * result;
}

constexpr std::expected<double, ArithmeticError> sqrt_approx(double x) {
  if (x < 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  if (x == 0)
    return 0.0;
  double guess = x;
  for (int i = 0; i < 20; ++i) {
    guess = 0.5 * (guess + x / guess);
  }
  return guess;
}

constexpr std::expected<double, ArithmeticError> factorial_approx(double x) {
  if (x < 0 || x != static_cast<int>(x)) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  if (x > 20) {
    return std::unexpected(ArithmeticError::Overflow);
  }
  return factorial(static_cast<int>(x));
}

constexpr std::expected<double, ArithmeticError> pow_approx(double a,
                                                            double b) {
  if (a == 0 && b <= 0) {
    return std::unexpected(ArithmeticError::DomainError);
  }
  if (b == static_cast<int>(b)) {
    double result = 1.0;
    int n = static_cast<int>(b);
    if (n < 0) {
      a = 1.0 / a;
      n = -n;
    }
    while (n > 0) {
      if (n % 2 == 1)
        result *= a;
      a *= a;
      n /= 2;
    }
    return result;
  }

  auto ln_a = ln_approx(a);
  if (!ln_a)
    return ln_a;
  double exp_arg = b * *ln_a;
  double result = 1.0;
  double term = 1.0;
  for (int n = 1; abs(term) > 1e-10 && n < 50; ++n) {
    result += term;
    term *= exp_arg / n;
  }
  return result;
}

constexpr std::expected<double, ArithmeticError> safe_add(double a, double b) {
  if (a > 0 && b > std::numeric_limits<double>::max() - a) {
    return std::unexpected(ArithmeticError::Overflow);
  }
  if (a < 0 && b < std::numeric_limits<double>::lowest() - a) {
    return std::unexpected(ArithmeticError::Underflow);
  }
  return a + b;
}

constexpr std::expected<double, ArithmeticError> safe_sub(double a, double b) {
  if (b > 0 && a < std::numeric_limits<double>::lowest() + b) {
    return std::unexpected(ArithmeticError::Underflow);
  }
  if (b < 0 && a > std::numeric_limits<double>::max() + b) {
    return std::unexpected(ArithmeticError::Overflow);
  }
  return a - b;
}

constexpr std::expected<double, ArithmeticError> safe_mul(double a, double b) {
  if (abs(a) > 1 && abs(b) > std::numeric_limits<double>::max() / abs(a)) {
    return std::unexpected(ArithmeticError::Overflow);
  }
  return a * b;
}

constexpr std::expected<double, ArithmeticError> safe_div(double a, double b) {
  if (abs(b) < 1e-10) {
    return std::unexpected(ArithmeticError::DivisionByZero);
  }
  return a / b;
}

enum class OpType {
  Add,
  Sub,
  Mul,
  Div,
  Pow,
  Sin,
  Cos,
  Tan,
  Cot,
  Sec,
  Csc,
  Asin,
  Acos,
  Atan,
  Acot,
  Asec,
  Acsc,
  Log,
  Ln,
  Sqrt,
  Factorial,
  Negate,
  LeftParen
};

template <Arithmetic T> struct Operation {
  OpType type;
  constexpr std::expected<T, ArithmeticError> apply(T a, T b = 0) const {
    switch (type) {
    case OpType::Add:
      return safe_add(a, b);
    case OpType::Sub:
      return safe_sub(a, b);
    case OpType::Mul:
      return safe_mul(a, b);
    case OpType::Div:
      return safe_div(a, b);
    case OpType::Pow:
      return pow_approx(a, b);
    case OpType::Sin:
      return sin_approx(a);
    case OpType::Cos:
      return cos_approx(a);
    case OpType::Tan:
      return tan_approx(a);
    case OpType::Cot:
      return cot_approx(a);
    case OpType::Sec:
      return sec_approx(a);
    case OpType::Csc:
      return csc_approx(a);
    case OpType::Asin:
      return asin_approx(a);
    case OpType::Acos:
      return acos_approx(a);
    case OpType::Atan:
      return atan_approx(a);
    case OpType::Acot:
      return acot_approx(a);
    case OpType::Asec:
      return asec_approx(a);
    case OpType::Acsc:
      return acsc_approx(a);
    case OpType::Log:
      return log10_approx(a);
    case OpType::Ln:
      return ln_approx(a);
    case OpType::Sqrt:
      return sqrt_approx(a);
    case OpType::Factorial:
      return factorial_approx(a);
    case OpType::Negate:
      return -a;
    default:
      return std::unexpected(ArithmeticError::UnknownOperator);
    }
  }
  constexpr bool is_unary() const {
    return type >= OpType::Sin && type <= OpType::Negate;
  }
};

template <Arithmetic T>
constexpr std::array<std::pair<std::string_view, Operation<T>>, 22>
    operation_map = {{{"+", Operation<T>{OpType::Add}},
                      {"-", Operation<T>{OpType::Sub}},
                      {"*", Operation<T>{OpType::Mul}},
                      {"/", Operation<T>{OpType::Div}},
                      {"^", Operation<T>{OpType::Pow}},
                      {"sin", Operation<T>{OpType::Sin}},
                      {"cos", Operation<T>{OpType::Cos}},
                      {"tan", Operation<T>{OpType::Tan}},
                      {"cot", Operation<T>{OpType::Cot}},
                      {"sec", Operation<T>{OpType::Sec}},
                      {"csc", Operation<T>{OpType::Csc}},
                      {"asin", Operation<T>{OpType::Asin}},
                      {"acos", Operation<T>{OpType::Acos}},
                      {"atan", Operation<T>{OpType::Atan}},
                      {"acot", Operation<T>{OpType::Acot}},
                      {"asec", Operation<T>{OpType::Asec}},
                      {"acsc", Operation<T>{OpType::Acsc}},
                      {"log", Operation<T>{OpType::Log}},
                      {"ln", Operation<T>{OpType::Ln}},
                      {"sqrt", Operation<T>{OpType::Sqrt}},
                      {"!", Operation<T>{OpType::Factorial}},
                      {"~", Operation<T>{OpType::Negate}}}};

template <Arithmetic T>
constexpr std::array<std::pair<std::string_view, int>, 22> precedence = {
    {{"+", 2},    {"-", 2},    {"*", 3},    {"/", 3},    {"^", 4},
     {"sin", 5},  {"cos", 5},  {"tan", 5},  {"cot", 5},  {"sec", 5},
     {"csc", 5},  {"asin", 5}, {"acos", 5}, {"atan", 5}, {"acot", 5},
     {"asec", 5}, {"acsc", 5}, {"log", 5},  {"ln", 5},   {"sqrt", 5},
     {"!", 5},    {"~", 5}}};

template <Arithmetic T>
constexpr std::optional<Operation<T>> find_operation(std::string_view op) {
  for (const auto &[name, operation] : operation_map<T>) {
    if (name == op)
      return operation;
  }
  return std::nullopt;
}

template <Arithmetic T> constexpr bool is_left_associative(OpType type) {
  return type != OpType::Pow;
}

template <Arithmetic T> struct Token {
  enum class Type { Number, Operator, LeftParen, RightParen, Function };
  Type type;
  std::variant<T, Operation<T>, std::string_view> value;
};

template <Arithmetic T> class Parser {
public:
  static std::expected<std::vector<Token<T>>, ArithmeticError>
  parse(std::string_view expression) {
    std::vector<Token<T>> tokens;
    bool expect_operand = true;
    size_t i = 0;

    while (i < expression.size()) {
      if (expression[i] == ' ') {
        i++;
        continue;
      }

      if (expect_operand && (expression[i] == '+' || expression[i] == '-')) {
        if (expression[i] == '-') {
          tokens.push_back(
              Token<T>{Token<T>::Type::Operator, Operation<T>{OpType::Negate}});
        }
        i++;
        continue;
      }

      if (expression[i] == '(') {
        tokens.push_back(Token<T>{Token<T>::Type::LeftParen, {}});
        expect_operand = true;
        i++;
      } else if (expression[i] == ')') {
        tokens.push_back(Token<T>{Token<T>::Type::RightParen, {}});
        expect_operand = false;
        i++;
      } else {
        size_t op_parse_start = i;
        std::string op_str;

        while (i < expression.size() &&
               (std::isalpha(expression[i]) || expression[i] == '!' ||
                expression[i] == '~')) {
          op_str += expression[i];
          i++;
        }

        if (!op_str.empty()) {
          auto op = find_operation<T>(op_str);
          if (op) {
            tokens.push_back(Token<T>{Token<T>::Type::Operator, *op});
            expect_operand = true;
          } else {
            i = op_parse_start;
          }
        }

        if (i == op_parse_start) {
          if (!expect_operand &&
              (expression[i] == '+' || expression[i] == '-' ||
               expression[i] == '*' || expression[i] == '/' ||
               expression[i] == '^')) {
            std::string_view bin_op_sv(&expression[i], 1);
            auto current_op_opt = find_operation<T>(bin_op_sv);
            if (!current_op_opt) {
              return std::unexpected(ArithmeticError::UnknownOperator);
            }
            tokens.push_back(
                Token<T>{Token<T>::Type::Operator, *current_op_opt});
            expect_operand = true;
            i++;
          } else if (expect_operand) {
            size_t num_start = i;
            bool has_sign_component =
                (expression[num_start] == '+' || expression[num_start] == '-');

            size_t num_end = num_start;
            if (num_end < expression.size() &&
                (expression[num_end] == '+' || expression[num_end] == '-')) {
              num_end++;
            }
            bool found_digit = false;
            while (num_end < expression.size() &&
                   (std::isdigit(expression[num_end]) ||
                    expression[num_end] == '.' || expression[num_end] == ',')) {
              if (std::isdigit(expression[num_end]))
                found_digit = true;
              num_end++;
            }

            if (!found_digit && has_sign_component &&
                num_end == num_start + 1) {
              return std::unexpected(ArithmeticError::InvalidInput);
            }
            if (num_end == num_start) {
              return std::unexpected(ArithmeticError::InvalidInput);
            }

            auto num_str = expression.substr(num_start, num_end - num_start);
            auto num = calc::parse_number<T>(num_str);
            if (!num)
              return std::unexpected(num.error());

            tokens.push_back(Token<T>{Token<T>::Type::Number, num.value()});
            expect_operand = false;
            i = num_end;
          } else {
            return std::unexpected(ArithmeticError::InvalidExpression);
          }
        }
      }
    }

    return tokens;
  }
};

template <Arithmetic T> class Evaluator {
public:
  static std::expected<T, ArithmeticError>
  evaluate(const std::vector<Token<T>> &tokens) {
    std::stack<T> values;
    std::stack<Operation<T>> ops;

    for (const auto &token : tokens) {
      switch (token.type) {
      case Token<T>::Type::Number:
        values.push(std::get<T>(token.value));
        break;
      case Token<T>::Type::Operator: {
        const auto &op = std::get<Operation<T>>(token.value);
        if (op.is_unary()) {
          ops.push(op);
        } else {
          while (!ops.empty() && ops.top().type != OpType::LeftParen) {
            const auto &top_stack_op = ops.top();
            int top_stack_op_precedence =
                precedence<T>[static_cast<size_t>(top_stack_op.type)].second;
            bool top_stack_op_is_left_assoc =
                is_left_associative<T>(top_stack_op.type);
            int current_op_precedence =
                precedence<T>[static_cast<size_t>(op.type)].second;

            if ((top_stack_op_is_left_assoc &&
                 top_stack_op_precedence >= current_op_precedence) ||
                (!top_stack_op_is_left_assoc &&
                 top_stack_op_precedence > current_op_precedence)) {
              auto result = apply_operator(values, ops);
              if (!result)
                return result;
              values.push(result.value());
            } else {
              break;
            }
          }
          ops.push(op);
        }
        break;
      }
      case Token<T>::Type::LeftParen:
        ops.push({OpType::LeftParen});
        break;
      case Token<T>::Type::RightParen:
        while (!ops.empty() && ops.top().type != OpType::LeftParen) {
          auto result = apply_operator(values, ops);
          if (!result)
            return result;
          values.push(result.value());
        }
        if (ops.empty()) {
          return std::unexpected(ArithmeticError::UnbalancedParentheses);
        }
        ops.pop();
        break;
      default:
        return std::unexpected(ArithmeticError::InvalidExpression);
      }
    }

    while (!ops.empty()) {
      if (ops.top().type == OpType::LeftParen) {
        return std::unexpected(ArithmeticError::UnbalancedParentheses);
      }
      auto result = apply_operator(values, ops);
      if (!result)
        return result;
      values.push(result.value());
    }

    if (values.size() != 1) {
      return std::unexpected(ArithmeticError::InvalidExpression);
    }

    return values.top();
  }

private:
  static std::expected<T, ArithmeticError>
  apply_operator(std::stack<T> &values, std::stack<Operation<T>> &ops) {
    if (ops.empty()) {
      return std::unexpected(ArithmeticError::InvalidExpression);
    }

    auto op = ops.top();
    ops.pop();

    if (op.is_unary()) {
      if (values.empty()) {
        return std::unexpected(ArithmeticError::InvalidExpression);
      }
      T a = values.top();
      values.pop();
      return op.apply(a);
    } else {
      if (values.size() < 2) {
        return std::unexpected(ArithmeticError::InvalidExpression);
      }
      T b = values.top();
      values.pop();
      T a = values.top();
      values.pop();
      return op.apply(a, b);
    }
  }
};

} // namespace detail

template <Arithmetic T>
std::expected<T, ArithmeticError> calculate(std::string_view expression) {
  try {
    auto tokens = detail::Parser<T>::parse(expression);
    if (!tokens)
      return std::unexpected(tokens.error());
    return detail::Evaluator<T>::evaluate(*tokens);
  } catch (...) {
    return std::unexpected(ArithmeticError::InvalidExpression);
  }
}

} // namespace calc