#include "seal/seal.h"
#include <iostream>
#include <vector>
#include <iomanip> // chronoライブラリのインクルード

using namespace std;
using namespace seal;

// 暗号化された行列の掛け算
vector<Ciphertext> encrypted_matrix_multiplication_batch(
    const vector<Ciphertext> &matrix1, const vector<Ciphertext> &matrix2,
    Evaluator &evaluator, RelinKeys &relin_keys, size_t n)
{
    vector<Ciphertext> result(n * n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // 初期化：最初の部分積
            size_t result_index = i * n + j;
            evaluator.multiply(matrix1[i * n + 0], matrix2[0 * n + j], result[result_index]);
            
            for (size_t k = 1; k < n; ++k) {
                Ciphertext temp;
                evaluator.multiply(matrix1[i * n + k], matrix2[k * n + j], temp);
                evaluator.add_inplace(result[result_index], temp); // 部分積を加算
            }

            // 再線形化
            evaluator.relinearize_inplace(result[result_index], relin_keys);
        }
    }
    return result;
}

// 結果を復号化して表示する関数
void display_matrix_result_batch(const vector<Ciphertext> &encrypted_matrix,
                                 Decryptor &decryptor, BatchEncoder &batch_encoder, size_t n) {
    vector<uint64_t> decoded_result(batch_encoder.slot_count());

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Plaintext plain_result;
            decryptor.decrypt(encrypted_matrix[i * n + j], plain_result);
            batch_encoder.decode(plain_result, decoded_result);
            cout << decoded_result[0] << " "; // スロット0の値を表示
        }
        cout << endl;
    }
}

// 平文行列の初期化を動的に行う
vector<uint64_t> generate_plain_matrix(size_t n, bool random = false) {
    vector<uint64_t> matrix(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix[i * n + j] = random ? rand() % 10 + 1 : (i * n + j + 1); // 1から10のランダム値またはデフォルト値
        }
    }
    return matrix;
}

int main() {
    // SEALパラメータの設定
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 16384;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20)); // ビット設定

    SEALContext context(parms);
    
    // キーの生成
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    // 行列のサイズ（n×n）: 動的に設定可能
    size_t n = 16; // 行列サイズ

    //行列文字（trueならランダム、falseならシーケンス（１、２、３...））
    vector<uint64_t> matrix1_plain = generate_plain_matrix(256, false);
    vector<uint64_t> matrix2_plain = generate_plain_matrix(256, false);

    
    // *** 実行開始からストラッセン法終了までの計測開始 ***
    auto start_time = chrono::high_resolution_clock::now();
    
    // 平文行列を暗号化
    vector<Ciphertext> matrix1_encrypted(n * n);
    vector<Ciphertext> matrix2_encrypted(n * n);

    for (size_t i = 0; i < n * n; ++i) {
        Plaintext plain1, plain2;
        batch_encoder.encode(vector<uint64_t>(batch_encoder.slot_count(), matrix1_plain[i]), plain1);
        batch_encoder.encode(vector<uint64_t>(batch_encoder.slot_count(), matrix2_plain[i]), plain2);

        encryptor.encrypt(plain1, matrix1_encrypted[i]);
        encryptor.encrypt(plain2, matrix2_encrypted[i]);
    }

    // 暗号化された行列の掛け算
    vector<Ciphertext> encrypted_result = encrypted_matrix_multiplication_batch(
        matrix1_encrypted, matrix2_encrypted, evaluator, relin_keys, n);
    
    // *** 計測終了 ***
    auto end_time = chrono::high_resolution_clock::now();
    
    // 実行時間の表示
    cout << fixed << setprecision(3)
         << "Execution Time (Start to End of Matrix Calculation): "
         << chrono::duration<double>(end_time - start_time).count() << " s" << endl;

    return 0;
}
